/**
 * @file fdmUtils.h
 * @author Special (special.xuan@outlook.com)
 * @brief
 * @version 1.3.2
 * @date 2024-06-25
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <math.h>
#include <time.h>
#include <udf.h> // note: move udf.h to the last

// debug
// #define DEBUG_FDM

// node match tolerance
#define NODE_MATCH_TOLERANCE 4e-6

// allocate new memories and initialize
#define NEW_MEMORIES(n, t, s) n = (t *)malloc(s * sizeof(t)); memset(n, 0, s * sizeof(t))

static int iTime = 0;    // record the number of the time steps
static int idFSI = 0;   // record the id of the fsi faces, shown in fluent
static int idFluid = 0; // record the id of the fluid cell zone, shown in fluent

static int nNodeFluid = 0;              // number of nodes in fluid
static int nModeFluid = 0;              // number of modes in fluid
static int nModeFewer = 0;              // number of modes in fluid
static int nNodeStruct = 0;             // number of nodes in fluid
static int nModeStruct = 0;             // number of modes in fluid
static double *structModeDisp = NULL;   // structure modal displacement
static double *structStiff = NULL;      // structure stiffness
static double *structDamp = NULL;       // structure damping
static double *structLoad = NULL;       // structure forces
static real *structModeForce = NULL;    // modal force
static real *modeFreq = NULL;           // frequencies of the structure
static real *modeForce_ThisTime = NULL; // modal force
static real *modeForce_LastTime = NULL; // modal force
static real *modeDisp_ThisTime = NULL;  // modal-displacement,modal-velocity,acceleration repeat
static real *modeDisp_LastTime = NULL;  // modal-displacement,modal-velocity,acceleration repeat
static real *initVelocity = NULL;       // initial modal velocity

/**
 * @brief compare nodes in order of x, y, z
 *
 * @param a pointer of first node
 * @param b pointer of second node
 * @return int
 */
int cmp_node(const void *const a, const void *const b)
{
    const double *da = (double *)a, *db = (double *)b;
    int xCmp = *(da + 0) - *(db + 0) < -NODE_MATCH_TOLERANCE ? -1 : *(da + 0) - *(db + 0) > NODE_MATCH_TOLERANCE ? 1 : 0; // compare x
    int yCmp = *(da + 1) - *(db + 1) < -NODE_MATCH_TOLERANCE ? -1 : *(da + 1) - *(db + 1) > NODE_MATCH_TOLERANCE ? 1 : 0; // compare y
    int zCmp = *(da + 2) - *(db + 2) < -NODE_MATCH_TOLERANCE ? -1 : *(da + 2) - *(db + 2) > NODE_MATCH_TOLERANCE ? 1 : 0; // compare z
    return xCmp ? xCmp : yCmp ? yCmp : zCmp ? zCmp : 0; // if x, then y, then z, then equal
}

/**
 * @brief sequencial search
 *
 * @param key key element
 * @param base base elements
 * @param nElements number of elements
 * @param sElemets size of elements
 * @return double*
 */
double *ssearch(const double *const key, const double *const base, const int nElements, const int sElements)
{
    double D = 0, *pResult = NULL;
    for (int i = 0; i < nElements; i++)
    {
        D = 0;
        for (int j = 0; j < ND_ND; j++)
            D += (key[j] - base[i * sElements + j]) * (key[j] - base[i * sElements + j]);
        D = sqrt(D);

        if (D < NODE_MATCH_TOLERANCE)
        {
            pResult = (double *)base + i * sElements;
            break;
        }
    }
    return pResult;
}

/**
 * @brief free all memories
 * 
 * @return int 
 */
void free_fdm_memories()
{
    free(modeFreq);
    free(modeForce_ThisTime);
    free(modeForce_LastTime);
    free(structModeForce);
    free(modeDisp_ThisTime);
    free(modeDisp_LastTime);
    free(initVelocity);

    Message0("UDF[Host]: All memories cleared!\n");
}

#ifdef DEBUG_FDM

/**
 * @brief output debug information
 * 
 * @param debugMessage debug information
 * @return int 
 */
int debug_file_print(const char *const debugMessage)
{
    static int debugFilePrintCount = 0;
    char debugFileName[30] = {0}; // output file name
    if I_AM_NODE_HOST_P
        sprintf(debugFileName, "debug/debugHost.csv");
    else
        sprintf(debugFileName, "debug/debugNode%d.csv", myid); // output file name of each node process
    FILE *fpDebug = fopen(debugFileName, "a");                 // open output file in write mode

    if (debugFilePrintCount == 0)
    {
        time_t timep;
        time(&timep);
        fprintf(fpDebug, "Calculation starts with %d processors, at %s", compute_node_count, asctime(localtime(&timep)));
    }
    fprintf(fpDebug, "%s", debugMessage);
    fclose(fpDebug);
    debugFilePrintCount++;
    return 0;
}
#endif

#if !RP_NODE
/**
 * @brief input number of nodes and modes of structure
 *
 * @return int
 */
int read_struct_nodes_modes()
{
    double nNodeFromFile = 0, nModeFromFile = 0;     // number of nodes, number of modes
    char inFileName[30] = "mode/StructNodeCoor.csv"; // input file name
    FILE *fpInput = fopen(inFileName, "r");          // open file read only
    if (fpInput == NULL)                             // file not exists print error
    {
        Message("UDF[Host]: Structure: Error: No Struct files.\n");
        return 1;
    }
    if (fscanf(fpInput, "%lf,", &nNodeFromFile) <= 0 || // ignore first number
        fscanf(fpInput, "%lf,", &nNodeFromFile) <= 0 || // input number of nodes
        fscanf(fpInput, "%lf,", &nModeFromFile) <= 0)   // input number of modes
    {
        Message("UDF[Host]: Structure: Error: Lack of Struct variables.\n");
        fclose(fpInput);
        return 2;
    }
    nNodeStruct = (int)nNodeFromFile; // number of nodes of structure
    nModeStruct = (int)nModeFromFile; // number of modes of structure
    fclose(fpInput);

    return 0;
}

/**
 * @brief input node coordinates and modal displacement
 *
 * @param structModeDisp node coordinates and modal displacements
 * @param row row of matrix
 * @param column column of matrix
 * @param iMode which mode to input
 * @return int
 */
int read_struct_coor_mode(double *const structModeDisp,
                          const int row,
                          const int column,
                          const int iMode)
{
    char inFileName[30] = {0}, line_buf[256] = {0}; // input file name, ignore first line
    if (iMode == 0)
        sprintf(inFileName, "mode/StructNodeCoor.csv"); // read coordinate at iMode == 0
    else
        sprintf(inFileName, "mode/StructNodeDisp%d.csv", iMode); // read modal displacement

    FILE *fpInput = fopen(inFileName, "r"); // open file read only
    if (fpInput == NULL)                    // file not exists print error
    {
        Message("UDF[Host]: Error: No Struct files.\n");
        return 1;
    }

    fgets(line_buf, 256, fpInput); // ignore first line
    for (int i = 0; i < row; i++)  // fill the array of node coordinates and modal displacements
        if (fscanf(fpInput, "%lf,", structModeDisp + i * column + (iMode - 1) * ND_ND + 0) <= 0 ||
            fscanf(fpInput, "%lf,", structModeDisp + i * column + (iMode - 1) * ND_ND + 1) <= 0 ||
            fscanf(fpInput, "%lf,", structModeDisp + i * column + (iMode - 1) * ND_ND + 2) <= 0) // data missing, print error
        {
            Message("UDF[Host]: Error: Lack of variables in Mode %d, Node %d.\n", iMode, i);
            fclose(fpInput);
            return 2;
        }

    fclose(fpInput);

    return 0;
}

/**
 * @brief input structure stiffness and damping matrices
 *
 * @param structStiff structure stiffness
 * @param structDamp structure damping
 * @param row row of matrices
 * @param column column of matrices
 * @return int
 */
int read_struct_stiff_damp_load(double *const structStiff,
                                double *const structDamp,
                                double *const structLoad,
                                const int row,
                                const int column)
{
    char line_buf[256] = {0};                           // ignore first line
    FILE *fpInput = fopen("mode/StructStiff.csv", "r"); // open file read only
    if (fpInput == NULL)                                // file not exists print error
    {
        Message("UDF[Host]: Error: No StructStiff file.\n");
        return 1;
    }
    fgets(line_buf, 256, fpInput); // ignore first line
    for (int i = 0; i < row; i++)  // fill the array of node coordinates and modal displacements
        if (fscanf(fpInput, "%lf,", structStiff + i * column + 0) <= 0 ||
            fscanf(fpInput, "%lf,", structStiff + i * column + 1) <= 0 ||
            fscanf(fpInput, "%lf,", structStiff + i * column + 2) <= 0) // data missing, print error
        {
            Message("UDF[Host]: Error: Lack of stiffness in Node %d.\n", i);
            fclose(fpInput);
            return 2;
        }
    fclose(fpInput);

    fpInput = fopen("mode/StructDamp.csv", "r"); // open file read only
    if (fpInput == NULL)                         // file not exists print error
    {
        Message("UDF[Host]: Error: No StructDamp file.\n");
        return 1;
    }
    fgets(line_buf, 256, fpInput); // ignore first line
    for (int i = 0; i < row; i++)  // fill the array of node coordinates and modal displacements
        if (fscanf(fpInput, "%lf,", structDamp + i * column + 0) <= 0 ||
            fscanf(fpInput, "%lf,", structDamp + i * column + 1) <= 0 ||
            fscanf(fpInput, "%lf,", structDamp + i * column + 2) <= 0) // data missing, print error
        {
            Message("UDF[Host]: Error: Lack of damp in Node %d.\n", i);
            fclose(fpInput);
            return 2;
        }
    fclose(fpInput);

    fpInput = fopen("mode/StructLoad.csv", "r"); // open file read only
    if (fpInput == NULL)                         // file not exists print error
    {
        Message("UDF[Host]: Error: No StructLoad file.\n");
        return 1;
    }
    fgets(line_buf, 256, fpInput); // ignore first line
    for (int i = 0; i < row; i++)  // fill the array of node coordinates and modal displacements
        if (fscanf(fpInput, "%lf,", structLoad + i * column + 0) <= 0 ||
            fscanf(fpInput, "%lf,", structLoad + i * column + 1) <= 0 ||
            fscanf(fpInput, "%lf,", structLoad + i * column + 2) <= 0) // data missing, print error
        {
            Message("UDF[Host]: Error: Lack of damp in Node %d.\n", i);
            fclose(fpInput);
            return 2;
        }
    fclose(fpInput);

    return 0;
}

/**
 * @brief input size of node coordinates and modal displacement
 *
 * @return int
 */
int read_fluid_nodes_modes()
{
    double nNodeFromFile = 0, nModeFromFile = 0;          // number of nodes, number of modes
    FILE *fpInput = fopen("mode/FluidNodeCoor.csv", "r"); // open the file in the read-only mode
    if (fpInput == NULL)                                  // file not exists, print error
    {
        Message("UDF[Host]: Error: No file.\n");
        return 1;
    }
    if (fscanf(fpInput, "%lf,", &nNodeFromFile) <= 0 || // ignore first number
        fscanf(fpInput, "%lf,", &nNodeFromFile) <= 0 || // input number of nodes
        fscanf(fpInput, "%lf,", &nModeFromFile) <= 0)   // input number of modes
    {
        Message("UDF[Host]: Error: Lack of variables.\n");
        fclose(fpInput);
        return 2;
    }
    nNodeFluid = (int)nNodeFromFile; // rows of node coordinates and modal displacement
    nModeFluid = (int)nModeFromFile; // columns of node coordinates and modal displacement
    fclose(fpInput);

    return 0;
}

/**
 * @brief input node coordinates and modal displacement
 *
 * @param nodeCoorDisp node coordinates and modal displacements
 * @param modeFreq mode frequency
 * @param iMode which mode to input
 * @return int
 */
int read_fluid_coor_mode(double *const nodeCoorDisp, real *const modeFreq, const int iMode)
{
    const int row = nNodeFluid, column = (nModeFluid + 1) * ND_ND;
    char inFileName[20] = {0}, line_buf[256] = {0}; // input file name, ignore first line
    if (iMode == 0)
        sprintf(inFileName, "mode/FluidNodeCoor.csv"); // read coordinate at first time
    else
        sprintf(inFileName, "mode/FluidNodeDisp%d.csv", iMode); // read modal displacement

    FILE *fpInput = fopen(inFileName, "r"); // open the file in the read-only mode
    if (fpInput == NULL)                    // file not exists, print error
    {
        Message("UDF[Host]: Error: No file.\n");
        return 1;
    }

    if (iMode > 0)
        if (fscanf(fpInput, "%lf,", modeFreq + iMode - 1) <=0) // read mode frequency
        {
            Message("UDF[Host]: Error: Lack of frequency in Mode %d.\n", iMode);
            fclose(fpInput);
            return 2;
        }

    fgets(line_buf, 256, fpInput); // ignore first line
    for (int i = 0; i < row; i++)  // fill the array of node coordinates and modal displacements
        if (fscanf(fpInput, "%lf,", nodeCoorDisp + i * column + iMode * ND_ND + 0) <= 0 ||
            fscanf(fpInput, "%lf,", nodeCoorDisp + i * column + iMode * ND_ND + 1) <= 0 ||
            fscanf(fpInput, "%lf,", nodeCoorDisp + i * column + iMode * ND_ND + 2) <= 0) // data missing, print error
        {
            Message("UDF[Host]: Error: Lack of variables in Mode %d, Node %d.\n", iMode, i);
            fclose(fpInput);
            return 2;
        }

    fclose(fpInput);
    return 0;
}

/**
 * @brief
 *
 * @return int
 */
int read_paramater()
{
    FILE *fpInput = fopen("mode/FluidPara.csv", "r"); // open the file in the read-only mode
    if (fpInput == NULL)                              // file not exists, print error
    {
        Message("UDF[Host]: Error: No file.\n");
        return 1;
    }
    if (fscanf(fpInput, "%d,\n", &idFSI) <= 0 ||
        fscanf(fpInput, "%d,\n", &idFluid) <= 0)
    {
        Message("UDF[Host]: Error: Lack of fluid id and fsi id.\n");
        fclose(fpInput);
        return 2;
    }
    for (int i = 0; i < nModeFluid; i++)
        if (fscanf(fpInput, "%lf,\n", initVelocity + i) <= 0) // data missing, print error
        {
            Message("UDF[Host]: Error: Lack of initial velocity in Mode %d.\n", i);
            fclose(fpInput);
            return 2;
        }
    Message("UDF[Host]: FSI id: %d\n           Fluid id : %d\n", idFSI, idFluid);

    fclose(fpInput);
    return 0;
}

/**
 * @brief Get structural modal force
 *
 * @return int
 */
int get_struct_mode_force()
{
    if (iTime <= 1) // skip first time step
        return 0;

    const real Dis = modeDisp_LastTime[nModeFluid * 3 + 0];         // modal displacement last time
    const real Vel = modeDisp_LastTime[nModeFluid * 3 + 1];         // modal velocity last time
    double *Stiff = (double *)malloc(nModeStruct * sizeof(double)); // Stiffness Force
    memset(Stiff, 0, nModeStruct * sizeof(double));
    double *Damp = (double *)malloc(nModeStruct * sizeof(double)); // Damp Force
    memset(Damp, 0, nModeStruct * sizeof(double));
    double *Load = (double *)malloc(nModeStruct * sizeof(double)); // Load Force
    memset(Load, 0, nModeStruct * sizeof(double));

    for (int iMode = 0; iMode < nModeStruct; iMode++)
    {
        for (int iNode = 0; iNode < nNodeStruct; iNode++)
        {
            int M = iNode * nModeStruct * ND_ND + iMode * ND_ND, N = iNode * ND_ND; // indices
            // if (iMode == nModeStruct -1 && iNode == nNodeStruct - 1)
            //     Message("M = %d, N = %d\n", M, N);

            Stiff[iMode] += structStiff[N + 0] * structModeDisp[N + 0] * structModeDisp[M + 0] +
                            structStiff[N + 1] * structModeDisp[N + 1] * structModeDisp[M + 1] +
                            structStiff[N + 2] * structModeDisp[N + 2] * structModeDisp[M + 2]; // stiffness force = K * x * x
            Damp[iMode] += structDamp[N + 0] * structModeDisp[N + 0] * structModeDisp[M + 0] +
                           structDamp[N + 1] * structModeDisp[N + 1] * structModeDisp[M + 1] +
                           structDamp[N + 2] * structModeDisp[N + 2] * structModeDisp[M + 2]; // damp force = C * x * x
            Load[iMode] += structLoad[N + 0] * structModeDisp[M + 0] +
                           structLoad[N + 1] * structModeDisp[M + 1] +
                           structLoad[N + 2] * structModeDisp[M + 2]; // load force = F * x
        }
        structModeForce[iMode] = (real)(Dis * Stiff[iMode] + Vel * Damp[iMode] + Load[iMode]); // f = F * x + D * K * x * x + V * C * x * x
    }

    free(Stiff); // clear memories
    free(Damp);
    free(Load);
    return 0;
}

/**
 * @brief wilson theta method
 *
 * @param dTime time gap
 * @return int
 */
int wilson_theta(const real dTime)
{
    static const real Mass = 1;    // modal mass = 1, which means normalization of principal mass
    static const real Damp = 0;    // modal damping
    static const real Theta = 1.4; // wilson-theta method, the value of the theta

    if (iTime <= 1) // initiate wilson-theta method
    {
        for (int i = 0; i < nModeFluid; i++)
        {
            modeDisp_ThisTime[i * 3 + 0] = 0;                    // D0 = 0
            modeDisp_ThisTime[i * 3 + 1] = initVelocity[i];      // V0
            modeDisp_ThisTime[i * 3 + 2] = modeForce_ThisTime[i]; // A0 = F0
        }
    }
    else
    {
        for (int i = 0; i < nModeFluid; i++)
        {
            real dis = 0, vel = 0, acc = 0;
            const real disLast = modeDisp_LastTime[i * 3 + 0], // displacement last time
                velLast = modeDisp_LastTime[i * 3 + 1],        // velocity last time
                accLast = modeDisp_LastTime[i * 3 + 2];        // acceleration last time

            // calculate modal displacement, acceleration, velocity using wilson-theta method
            dis = modeForce_LastTime[i] + Theta * (modeForce_ThisTime[i] - modeForce_LastTime[i]);
            dis += Mass * (6 / SQR(Theta * dTime) * disLast + 6 / (Theta * dTime) * velLast + 2 * accLast);
            dis += Damp * (3 / (Theta * dTime) * disLast + 2 * velLast + 0.5 * Theta * dTime * accLast);
            dis /= 6 * Mass / SQR(Theta * dTime) + 3 * Damp / (Theta * dTime) + SQR(2 * M_PI * modeFreq[i]);
            // calculate acceleration
            acc = 6 / (SQR(Theta * dTime) * Theta) * (dis - disLast) - 6 / (SQR(Theta) * dTime) * velLast + (1 - 3 / Theta) * accLast;
            // calculate velocity
            vel = velLast + 0.5 * dTime * (acc + accLast);
            // calculate displacement
            dis = disLast + dTime * velLast + SQR(dTime) / 6 * (acc + 2 * accLast);

            modeDisp_ThisTime[i * 3 + 0] = dis; // displacement this time
            modeDisp_ThisTime[i * 3 + 1] = vel; // velocity this time
            modeDisp_ThisTime[i * 3 + 2] = acc; // displacement this time
        }
    }

    return 0;
}
#endif

#if !RP_HOST
/**
 * @brief fill node coordinate and modal displacement into UDMI
 *
 * @param nodeCoorDisp node coordinate and modal displacement
 * @return int
 */
int fill_modal_disp(const double *const nodeCoorDisp)
{
    Domain *pDomain; // pointer of domain
    cell_t pCell;    // pointer of cell
    Thread *pThread; // pointer of thread
    Node *pNode;     // pointer of node

    const int row = nNodeFluid, column = (nModeFluid + 1) * ND_ND;
    int iNode = 0, nodeCount = 0, nodeErrorCount = 0, UDMIColumn = column - ND_ND; // node index, total number of node, columns of UDMI
    double node_coor[ND_ND] = {0}, *thisNode = NULL;           // buffer of node coordinate, pointer of this node

    char outFileName[30] = {0};                          // output file name
    sprintf(outFileName, "FluidOutputNode%d.csv", myid); // output file name of each node process
    FILE *fpOutput = fopen(outFileName, "w+");           // open output file in write mode

    pDomain = Get_Domain(1);        // for single-phase flows, domain_id is 1 and Get_Domain(1) returns the fluid domain pointer
    thread_loop_c(pThread, pDomain) // loop thread in domain
    {
        begin_c_loop_int_ext(pCell, pThread) // loop cell in thread
        {
            c_node_loop(pCell, pThread, iNode) // loop node in cell
            {
                pNode = C_NODE(pCell, pThread, iNode); // get node pointer
                N_UDMI(pNode, UDMIColumn) = 0;         // last column of UDMI is 0 indicates not set, 1 indicates set
            }
        }
    end_c_loop_int_ext(pCell, pThread) // finish cell loop
    }
    thread_loop_c(pThread, pDomain) // loop thread in domain
    {
        begin_c_loop_int_ext(pCell, pThread) // loop cell in thread
        {
            c_node_loop(pCell, pThread, iNode) // loop node in cell
            {
                pNode = C_NODE(pCell, pThread, iNode); // get node pointer
                if (N_UDMI(pNode, UDMIColumn) != 0)    // if node set
                    continue;
                
                node_coor[0] = NODE_X(pNode); // get coordinate of this node
                node_coor[1] = NODE_Y(pNode);
                node_coor[2] = NODE_Z(pNode);

                thisNode = (double *)bsearch(node_coor, nodeCoorDisp, row, column * sizeof(double), cmp_node); // search for this node
                if (thisNode != NULL)                    // if this node found
                {
                    for (int i = 0; i < UDMIColumn; i++) // fill UDMI with modal displacement
                        N_UDMI(pNode, i) = thisNode[i + ND_ND];
                    N_UDMI(pNode, UDMIColumn) = 1; // this node set
                    nodeCount++;                   // count how many nodes in this process
                    continue;
                }

                thisNode = ssearch(node_coor, nodeCoorDisp, row, column);
                if (thisNode != NULL)                    // if this node found
                {
                    for (int i = 0; i < UDMIColumn; i++) // fill UDMI with modal displacement
                        N_UDMI(pNode, i) = thisNode[i + ND_ND];
                    fprintf(fpOutput, "node at (%20.10e, %20.10e, %20.10e), (%20.10e, %20.10e, %20.10e)\n", node_coor[0], node_coor[1], node_coor[2], thisNode[0], thisNode[1], thisNode[2]); // if this node not found, output error in file
                }
                else
                {
                    fprintf(fpOutput, "Error: no node at (%20.10e, %20.10e, %20.10e)\n", node_coor[0], node_coor[1], node_coor[2]); // if this node not found, output error in file
                    nodeErrorCount++;
                }
                N_UDMI(pNode, UDMIColumn) = 1; // this node set
                nodeCount++;                   // count how many nodes in this process
            }
        }
        end_c_loop_int_ext(pCell, pThread) // finish cell loop
        }
    fclose(fpOutput);                                                                // close file
    Message("UDF[Node]: Total Number of nodes in NODE %d is %d, Error nodes is %d\n", myid, nodeCount, nodeErrorCount); // print how many nodes in this process

    return 0;
}

/**
 * @brief get mode force object
 *
 * @param pDomain pointer of domain
 * @return int
 */
int get_fluid_mode_force(Domain *const pDomain)
{
    Thread *pThread = Lookup_Thread(pDomain, idFSI); // pointer of thread
    face_t pFace;                                    // pointer of face
    Node *pNode;                                     // pointer of node
    real NV_VEC(AreaVector) = {0},                   // normal vector of area
        NV_VEC(ForceVector) = {0},                   // pressure vector
        NV_VEC(ShearVector) = {0},                   // Shear vector
        Press = 0;                                   // pressure
    int nNodePerFace = 0, iNode = 0;                 // number of node per face, node index
    static int modeForceCount = 0;                   // modal force counter

#ifdef DEBUG_FDM
    real NV_VEC(AreaCentroid) = {0};
    char debugMessage[200] = {0};
#endif

    begin_f_loop(pFace, pThread) // begin face loop
    if (PRINCIPAL_FACE_P(pFace, pThread))
    {
        nNodePerFace = F_NNODES(pFace, pThread); // get number of nodes per face
        F_AREA(AreaVector, pFace, pThread);      // get normal vector
        Press = F_P(pFace, pThread);             // get pressure
        NV_VS(ShearVector, =, F_STORAGE_R_NV(pFace, pThread, SV_WALL_SHEAR), *, -1.0); // get shear vector
        NV_V_VS(ForceVector, =, ShearVector, +, AreaVector, *, Press); // calculate pressure vector
#ifdef DEBUG_FDM
        if (iTime == 2)
        {
            F_Centroid(AreaCentroid, pFace, pThread);
            sprintf(debugMessage, "%5d, %20.10e, %20.10e, %5d,\n",
                    pFace, ForceVector[0], ShearVector[0], nNodePerFace);
            debug_file_print(debugMessage);
        }
#endif
        f_node_loop(pFace, pThread, iNode) // begin node loop
        {
            pNode = F_NODE(pFace, pThread, iNode); // get this node
            for (int i = 0; i < nModeFluid; i++)   // for each mode
                for (int j = 0; j < ND_ND; j++)
                    modeForce_ThisTime[i] += 1 / (double)nNodePerFace * ForceVector[j] * N_UDMI(pNode, i * ND_ND + j);

            modeForceCount++; // modal force counter +1
        }
    }
    end_f_loop(pFace, pThread); // end face loop

    return 0;
}

/**
 * @brief move grid
 *
 * @param pDomain
 * @return int
 */
int move_grid(Domain *const pDomain)
{
    if (iTime <= 1) // skip first time step
        return 0;

    Thread *pThread = Lookup_Thread(pDomain, idFluid); // pointer of thread
    cell_t pCell;                                      // pointer of cell
    Node *pNode;                                       // pointer of node
    real NV_VEC(dispUpdate) = {0}, deltaDisp = 0;      // update displacement, displacement gap
    int iNode = 0;                                     // node index

    begin_c_loop_int_ext(pCell, pThread) // begin all cell loop
    {
        c_node_loop(pCell, pThread, iNode) // begin node loop
        {
            pNode = C_NODE(pCell, pThread, iNode); // get node pointer
            if (NODE_POS_NEED_UPDATE(pNode))       // if node not yet updated
            {
                memset(dispUpdate, 0, ND_ND * sizeof(real)); // clear update displacement
                for (int i = 0; i < nModeFluid; i++)
                {
                    deltaDisp = modeDisp_ThisTime[i * 3] - modeDisp_LastTime[i * 3]; // displacement gap
                    dispUpdate[0] += deltaDisp * N_UDMI(pNode, i * ND_ND + 0);       // update displacement is summation of displacement gap times each modal displacement
                    dispUpdate[1] += deltaDisp * N_UDMI(pNode, i * ND_ND + 1);
                    dispUpdate[2] += deltaDisp * N_UDMI(pNode, i * ND_ND + 2);
                }
                NV_V(NODE_COORD(pNode), +=, dispUpdate); // update node
                NODE_POS_UPDATED(pNode);                 // indicate node updated
            }
        }
    }
    end_c_loop_int_ext(pCell, pThread) // end cell loop

    return 0;
}
#endif
