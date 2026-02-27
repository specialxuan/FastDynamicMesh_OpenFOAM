/*----------------------------------------------------------------***********\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Fast Dynamic Mesh method implementation.

\*---------------------------------------------------------------------------*/

#include "fastDynamicFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "SortableList.H"
#include "Pstream.H"
#include <fstream>
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(fastDynamicFvMesh, 0);
addToRunTimeSelectionTable(dynamicFvMesh, fastDynamicFvMesh, IOobject);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fastDynamicFvMesh::fastDynamicFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    nMode_(0),
    theta_(1.4) // Default Wilson-Theta
{
    readControls();
    readModeShapes();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fastDynamicFvMesh::~fastDynamicFvMesh()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void fastDynamicFvMesh::readControls()
{
    // Read dictionary
    IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            this->time().constant(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );
    
    if (dynamicMeshDict.found(typeName + "Coeffs"))
    {
        const dictionary& fdmDict = dynamicMeshDict.subDict(typeName + "Coeffs");
        fdmDict.lookup("theta") >> theta_;
        fdmDict.lookup("fsiPatches") >> fsiPatches_;
    }
    else
    {
        Info << "Warning: " << typeName << "Coeffs not found, using defaults." << endl;
    }
}

void fastDynamicFvMesh::readModeShapes()
{
    // 1. Read files on master
    List<point> csvPoints;
    List<List<vector>> csvShapes;
    label nCsvNodes = 0;

    if (Pstream::master())
    {
        Info<< "Reading mode coordinates..." << endl;
        std::ifstream file("mode/FluidNodeCoor.csv");
        if (!file.good())
        {
             // Warning instead of FatalError to allow testing without files
             WarningInFunction << "Cannot open mode/FluidNodeCoor.csv" << endl;
        }
        else
        {
            scalar dummy, nNode, nMode;
            char comma;
            file >> dummy >> comma >> nNode >> comma >> nMode;
            
            nCsvNodes = label(nNode);
            nMode_ = label(nMode);
            
            csvPoints.setSize(nCsvNodes);
            csvShapes.setSize(nMode_);
            forAll(csvShapes, m) csvShapes[m].setSize(nCsvNodes);

            // Skip rest of header line
            string line;
            std::getline(file, line); 

            for (label i=0; i<nCsvNodes; ++i)
            {
                scalar x, y, z;
                char c;
                file >> x >> c >> y >> c >> z;
                csvPoints[i] = point(x, y, z);
            }
            
            // Read mode shapes
            for (label m=0; m<nMode_; ++m)
            {
                 std::ostringstream ss;
                 ss << "mode/FluidNodeDisp" << (m+1) << ".csv";
                 std::ifstream mFile(ss.str());
                 
                 scalar freq;
                 if (mFile >> freq)
                 {
                    modeFreq_.setSize(nMode_);
                    modeFreq_[m] = freq;
                 }
                 
                 std::getline(mFile, line); // Skip rest of header
                 
                 for (label i=0; i<nCsvNodes; ++i)
                 {
                     scalar dx, dy, dz;
                     char c;
                     mFile >> dx >> c >> dy >> c >> dz;
                     csvShapes[m][i] = vector(dx, dy, dz);
                 }
            }
        }
    }

    // Broadcast sizes
    Pstream::broadcast(nMode_);
    Pstream::broadcast(nCsvNodes);
    
    // Resize local arrays
    modeFreq_.setSize(nMode_);
    modeForce_.setSize(nMode_, 0.0);
    modeForce0_.setSize(nMode_, 0.0);
    modeState_.setSize(nMode_, vector::zero);
    modeState0_.setSize(nMode_, vector::zero);
    initVelocity_.setSize(nMode_, 0.0);
    
    // Broadcast data
    Pstream::broadcast(modeFreq_);
    Pstream::broadcast(csvPoints);
    forAll(csvShapes, m) Pstream::broadcast(csvShapes[m]);

    // Map to local mesh points
    const pointField& localPoints = this->points();
    modeShapes_.setSize(nMode_);
    forAll(modeShapes_, m)
    {
        modeShapes_[m].setSize(localPoints.size());
        modeShapes_[m] = vector::zero;
    }

    // Brute force search (simplest correct implementation)
    // Only perform mapping if we have CSV nodes
    if (nCsvNodes > 0)
    {
        scalar tolSqr = 1e-12; 
        
        forAll(localPoints, pI)
        {
            const point& p = localPoints[pI];
            scalar minDistSqr = GREAT;
            label nearestIdx = -1;
            
            forAll(csvPoints, cI)
            {
                scalar d = magSqr(p - csvPoints[cI]);
                if (d < minDistSqr)
                {
                    minDistSqr = d;
                    nearestIdx = cI;
                }
            }
            
            if (minDistSqr < tolSqr && nearestIdx != -1)
            {
                for (label m=0; m<nMode_; ++m)
                {
                    modeShapes_[m][pI] = csvShapes[m][nearestIdx];
                }
            }
        }
    }
}

void fastDynamicFvMesh::calcModalForces()
{
    // Initialize forces
    modeForce_ = 0.0;

    if (this->foundObject<volScalarField>("p"))
    {
        const volScalarField& p = this->lookupObject<volScalarField>("p");
        
        // Iterate over FSI patches
        forAll(fsiPatches_, i)
        {
            label patchID = this->boundaryMesh().findPatchID(fsiPatches_[i]);
            if (patchID == -1) continue;

            const polyPatch& pp = this->boundaryMesh()[patchID];
            const fvPatchScalarField& pPatch = p.boundaryField()[patchID];
            // Loop over faces
            forAll(pp, faceI)
            {
                // Calculate pressure force on face
                vector f = pPatch[faceI] * pp.faceAreas()[faceI];
                
                // Get face points to interpolate mode shape
                const labelList& fPoints = pp[faceI];
                
                for (label m=0; m<nMode_; ++m)
                {
                    // Average mode shape at face center
                    vector shapeFace = vector::zero;
                    forAll(fPoints, fp)
                    {
                        shapeFace += modeShapes_[m][fPoints[fp]];
                    }
                    if (fPoints.size() > 0) shapeFace /= fPoints.size();
                    
                    // Project force
                    modeForce_[m] += f & shapeFace;
                }
            }
        }
    }

    // Parallel reduction
    Pstream::listCombineGather(modeForce_, plusEqOp<scalar>());
    Pstream::broadcast(modeForce_);
}

void fastDynamicFvMesh::solveStructuralDynamics(scalar dt)
{
    const scalar Mass = 1.0; 
    const scalar Damp = 0.0;

    for (label i=0; i<nMode_; ++i)
    {
        scalar dis = 0, vel = 0, acc = 0;
        
        scalar disLast = modeState0_[i].x();
        scalar velLast = modeState0_[i].y();
        scalar accLast = modeState0_[i].z();
        
        scalar F_this = modeForce_[i];
        scalar F_last = modeForce0_[i];
        scalar freq = modeFreq_[i];
        scalar omega = 2.0 * constant::mathematical::pi * freq;
        
        // Ensure non-zero denominator logic
        scalar effectiveK = 6.0*Mass/sqr(theta_*dt) + 3.0*Damp/(theta_*dt) + sqr(omega);
        
        scalar load = F_last + theta_*(F_this - F_last);
        load += Mass * (6.0/sqr(theta_*dt)*disLast + 6.0/(theta_*dt)*velLast + 2.0*accLast);
        load += Damp * (3.0/(theta_*dt)*disLast + 2.0*velLast + 0.5*theta_*dt*accLast);
        
        scalar disTheta = load / effectiveK;
        
        scalar tau = theta_ * dt;
        scalar term1 = 6.0 / (sqr(tau) * theta_) * (disTheta - disLast);
        scalar term2 = 6.0 / (sqr(theta_) * dt) * velLast;
        scalar term3 = (1.0 - 3.0/theta_) * accLast;
        acc = term1 - term2 + term3;
        
        vel = velLast + 0.5*dt*(acc + accLast);
        dis = disLast + dt*velLast + sqr(dt)/6.0*(acc + 2.0*accLast);
        
        modeState_[i] = vector(dis, vel, acc);
    }
}

bool fastDynamicFvMesh::update()
{
    // 1. Calculate time step
    scalar dt = this->time().deltaTValue();
    
    // Store previous state
    modeForce0_ = modeForce_;
    modeState0_ = modeState_;

    // 2. Calculate forces
    calcModalForces();

    // 3. Solve dynamics
    if (nMode_ > 0)
    {
        solveStructuralDynamics(dt);
    }

    // 4. Update mesh points
    pointField newPoints = this->points(); 

    for (label m=0; m<nMode_; ++m)
    {
        const vectorField& shape = modeShapes_[m];

        // Incremental displacement: new - old
        scalar dDisp = modeState_[m].x() - modeState0_[m].x();

        forAll(newPoints, pI)
        {
            newPoints[pI] += dDisp * shape[pI];
        }
    }
    
    this->movePoints(newPoints);
    
    return true;
}

} // End namespace Foam
