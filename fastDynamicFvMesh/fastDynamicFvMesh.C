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
#include "fvcGrad.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "fluidThermo.H"
#include "transportModel.H"
#include "syncTools.H"
#include "SortableList.H"
#include "Pstream.H"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>

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
    theta_(1.4), // Default Wilson-Theta
    mappingTolerance_(4e-6),
    startupStepCount_(0)
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
        fdmDict.readIfPresent("theta", theta_);
        fdmDict.lookup("fsiPatches") >> fsiPatches_;
        fdmDict.readIfPresent("mappingTolerance", mappingTolerance_);
    }
    else
    {
        Info << "Warning: " << typeName << "Coeffs not found, using defaults." << endl;
    }
}

void fastDynamicFvMesh::readLegacyParameters(const fileName& modeDir)
{
    if (!Pstream::master())
    {
        return;
    }

    fileName paraPath = modeDir/"FluidPara.csv";
    std::ifstream paraFile(paraPath);

    if (!paraFile.good())
    {
        Info<< "Legacy parameter file " << paraPath
            << " not found; using zero initial modal velocities." << endl;
        return;
    }

    auto readCsvScalar = [&](scalar& value) -> bool
    {
        string line;

        while (std::getline(paraFile, line))
        {
            std::replace(line.begin(), line.end(), ',', ' ');
            std::stringstream ss(line);

            if (ss >> value)
            {
                return true;
            }
        }

        return false;
    };

    scalar legacyFsiId = -1;
    scalar legacyFluidId = -1;

    if (!readCsvScalar(legacyFsiId) || !readCsvScalar(legacyFluidId))
    {
        WarningInFunction
            << "Cannot read legacy ids from " << paraPath
            << ". Initial modal velocities remain zero." << endl;

        return;
    }

    Info<< "Read legacy FluidPara.csv ids (fsi=" << legacyFsiId
        << ", fluid=" << legacyFluidId
        << "). OpenFOAM uses patch names from dynamicMeshDict instead."
        << endl;

    for (label modeI = 0; modeI < nMode_; ++modeI)
    {
        scalar value = 0.0;

        if (!readCsvScalar(value))
        {
            WarningInFunction
                << "Missing initial velocity for mode " << modeI
                << " in " << paraPath
                << ". Remaining initial velocities stay at zero." << endl;

            break;
        }

        initVelocity_[modeI] = value;
    }
}

void fastDynamicFvMesh::readModeShapes()
{
    // 1. Read files on master
    List<point> csvPoints;
    List<List<vector>> csvShapes;
    label nCsvNodes = 0;
    fileName modeDir = this->time().path()/"mode";

    if (Pstream::parRun() && !exists(modeDir))
    {
        modeDir = this->time().path().path()/"mode";
    }

    if (Pstream::master())
    {
        Info<< "Reading mode coordinates..." << endl;

        fileName coorPath = modeDir/"FluidNodeCoor.csv";
        Info<< "Trying to open: " << coorPath << endl;

        std::ifstream file(coorPath);
        if (!file.good())
        {
             // Warning instead of FatalError to allow testing without files
             WarningInFunction << "Cannot open " << coorPath << endl;
        }
        else
        {
            scalar dummy, nNode, nMode;
            // Add robust parsing for header
            // Read line, replace commas with spaces, read numbers
            string line;
            std::getline(file, line);
            std::replace(line.begin(), line.end(), ',', ' ');
            std::stringstream ss(line);
            ss >> dummy >> nNode >> nMode;
            
            nCsvNodes = label(nNode);
            nMode_ = label(nMode);
            
            Info<< "  Found " << nCsvNodes << " nodes and " << nMode_ << " modes in CSV." << endl;

            modeFreq_.setSize(nMode_); // Resize here on master
            csvPoints.setSize(nCsvNodes);
            csvShapes.setSize(nMode_);
            forAll(csvShapes, m) csvShapes[m].setSize(nCsvNodes);

            // Skip rest of first line (already read)
            // std::getline(file, line); 

            for (label i=0; i<nCsvNodes; ++i)
            {
                scalar x, y, z;
                // Robust parsing for coordinates
                // Read 3 numbers separated by commas
                // file >> x >> c >> y >> c >> z; 
                // This fails if there are spaces around comma.
                // Better: read line, replace commas, read numbers
                std::getline(file, line);
                std::replace(line.begin(), line.end(), ',', ' ');
                std::stringstream ss2(line);
                ss2 >> x >> y >> z;
                
                csvPoints[i] = point(x, y, z);
            }
            
            // Read mode shapes
            for (label m=0; m<nMode_; ++m)
            {
                 fileName shapePath = modeDir/("FluidNodeDisp" + std::to_string(m+1) + ".csv");
                 std::ifstream mFile(shapePath);
                 
                 string line;
                 std::getline(mFile, line);
                 // First line is frequency
                 std::stringstream ss(line);
                 scalar freq;
                 if (ss >> freq)
                 {
                    modeFreq_[m] = freq;
                    Info<< "  Mode " << m << " Freq: " << freq << endl;
                 }
                 else
                 {
                    Info<< "  Warning: Failed to read frequency for mode " << m << endl;
                 }
                 
                 // Second line is header, skip?
                 // Wait, original logic read `freq` then skipped line.
                 // The file format:
                 // Frequency
                 // Header
                 // Data
                 
                 std::getline(mFile, line); // Skip header line
                 
                 for (label i=0; i<nCsvNodes; ++i)
                 {
                     std::getline(mFile, line);
                     std::replace(line.begin(), line.end(), ',', ' ');
                     std::stringstream ss2(line);
                     
                     scalar dx, dy, dz;
                     ss2 >> dx >> dy >> dz;
                     csvShapes[m][i] = vector(dx, dy, dz);
                 }
            }
        }
    }

    // Broadcast sizes
    Pstream::broadcast(nMode_);
    Pstream::broadcast(nCsvNodes);
    
    // Resize local arrays
    if (!Pstream::master()) modeFreq_.setSize(nMode_); // Only slaves need resize now

    modeForce_.setSize(nMode_, 0.0);
    modeForce0_.setSize(nMode_, 0.0);
    modeState_.setSize(nMode_, vector::zero);
    modeState0_.setSize(nMode_, vector::zero);
    initVelocity_.setSize(nMode_, 0.0);

    readLegacyParameters(modeDir);

    // Broadcast data
    Pstream::broadcast(modeFreq_);
    Pstream::broadcast(csvPoints);
    Pstream::broadcast(initVelocity_);

    // Resize csvShapes on all procs to ensure loop consistency
    csvShapes.setSize(nMode_);
    
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
    label mappedCount = 0;
    if (nCsvNodes > 0)
    {
        scalar tolSqr = sqr(mappingTolerance_);

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
                mappedCount++;
                for (label m=0; m<nMode_; ++m)
                {
                    modeShapes_[m][pI] = csvShapes[m][nearestIdx];
                }
            }
        }
    }
    Info<< "Mapped " << mappedCount << " points out of " << localPoints.size() << " local mesh points." << endl;
}

tmp<scalarField> fastDynamicFvMesh::patchDensity
(
    const label patchi,
    const scalar defaultRho
) const
{
    if (this->foundObject<volScalarField>("rho"))
    {
        const auto& rhoField = this->lookupObject<volScalarField>("rho");
        return tmp<scalarField>
        (
            new scalarField(rhoField.boundaryField()[patchi])
        );
    }

    return tmp<scalarField>
    (
        new scalarField(this->boundaryMesh()[patchi].size(), defaultRho)
    );
}

tmp<symmTensorField> fastDynamicFvMesh::devRhoReff
(
    const tensorField& gradUp,
    const label patchi,
    const scalar defaultRho
) const
{
    typedef incompressible::turbulenceModel icoTurbModel;
    typedef compressible::turbulenceModel cmpTurbModel;

    if (this->foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const auto& turb =
            this->lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        tmp<scalarField> tRho = patchDensity(patchi, defaultRho);

        return tmp<symmTensorField>
        (
            new symmTensorField
            (
                -tRho()*turb.nuEff(patchi)*devTwoSymm(gradUp)
            )
        );
    }
    else if (this->foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const auto& turb =
            this->lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        return tmp<symmTensorField>
        (
            new symmTensorField
            (
                -turb.muEff(patchi)*devTwoSymm(gradUp)
            )
        );
    }
    else if (this->foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const auto& thermo = this->lookupObject<fluidThermo>(fluidThermo::dictName);

        return tmp<symmTensorField>
        (
            new symmTensorField
            (
                -thermo.mu(patchi)*devTwoSymm(gradUp)
            )
        );
    }
    else if (this->foundObject<transportModel>("transportProperties"))
    {
        const auto& laminarT =
            this->lookupObject<transportModel>("transportProperties");

        tmp<scalarField> tRho = patchDensity(patchi, defaultRho);

        return tmp<symmTensorField>
        (
            new symmTensorField
            (
                -tRho()*laminarT.nu(patchi)*devTwoSymm(gradUp)
            )
        );
    }
    else
    {
        IOdictionary transportProperties
        (
            IOobject
            (
                "transportProperties",
                this->time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        const dimensionedScalar nu("nu", dimViscosity, transportProperties);
        tmp<scalarField> tRho = patchDensity(patchi, defaultRho);

        return tmp<symmTensorField>
        (
            new symmTensorField
            (
                -tRho()*nu.value()*devTwoSymm(gradUp)
            )
        );
    }
}

void fastDynamicFvMesh::calcModalForces()
{
    // Initialize forces
    modeForce_ = 0.0;

    // Default density (kinematic = 1.0)
    // Hardcoded for water validation case
    scalar rho = 1000.0;

    // Check if we have rho field (compressible)
    if (this->foundObject<volScalarField>("rho"))
    {
        // For compressible, forces should naturally integrate p (which is P_abs).
        // BUT, OpenFOAM compressible p is usually Pa.
        // Incompressible p is m^2/s^2.
        // We assume here the user wants forces in Newtons.
        // If "rho" field exists, p is likely Pa, so we don't multiply by rho.
        // So rho factor remains 1.0.
        rho = 1.0;
    }
    else
    {
        // Incompressible: p is p_kinematic. Need rhoRef.
        IOdictionary transportProperties
        (
            IOobject
            (
                "transportProperties",
                this->time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        if (transportProperties.found("rho"))
        {
             transportProperties.lookup("rho") >> rho;
        }
    }

    if (Pstream::master())
    {
        Info<< "DEBUG: Using rho=" << rho << endl;
    }

    tmp<volTensorField> tGradU;
    const volTensorField* gradUPtr = nullptr;

    if (this->foundObject<volVectorField>("U"))
    {
        const volVectorField& U = this->lookupObject<volVectorField>("U");
        tGradU = fvc::grad(U);
        gradUPtr = &tGradU();
    }
    else
    {
        WarningInFunction
            << "Velocity field 'U' not found; wall shear contribution "
            << "to modal forces will be skipped." << endl;
    }

    if (this->foundObject<volScalarField>("p"))
    {
        const volScalarField& p = this->lookupObject<volScalarField>("p");

        // Debug p range
        scalar pMin = gMin(p);
        scalar pMax = gMax(p);
        if (Pstream::master())
        {
             Info << "DEBUG: p range [" << pMin << ", " << pMax << "]" << endl;
        }

        // Iterate over FSI patches
        forAll(fsiPatches_, i)
        {
            label patchID = this->boundaryMesh().findPatchID(fsiPatches_[i]);
            if (patchID == -1)
            {
                Info<< "Warning: FSI patch " << fsiPatches_[i] << " not found!" << endl;
                continue;
            }

            const polyPatch& pp = this->boundaryMesh()[patchID];
            const fvPatchScalarField& pPatch = p.boundaryField()[patchID];
            const vectorField faceAreas(pp.faceAreas());

            const symmTensorField* devStressPtr = nullptr;
            tmp<symmTensorField> tDevStress;

            if (gradUPtr)
            {
                const tensorField& gradPatch = gradUPtr->boundaryField()[patchID];
                tDevStress = devRhoReff(gradPatch, patchID, rho);
                devStressPtr = &tDevStress();
            }

            // Loop over faces
            forAll(pp, faceI)
            {
                const vector& areaVec = faceAreas[faceI];
                vector faceForce = rho * pPatch[faceI] * areaVec;

                if (devStressPtr)
                {
                    faceForce += areaVec & (*devStressPtr)[faceI];
                }

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
                    if (fPoints.size() > 0)
                    {
                        shapeFace /= fPoints.size();
                    }

                    // Project force
                    modeForce_[m] += faceForce & shapeFace;
                }
            }
        }
    }

    // Parallel reduction
    Pstream::listCombineGather(modeForce_, plusEqOp<scalar>());
    Pstream::broadcast(modeForce_);

    if (Pstream::master())
    {
         Info << "DEBUG: rho=" << rho << endl;
         Info << "DEBUG: Mode Forces (0-4): ";
         for(label i=0; i<min(5, nMode_); ++i) Info << modeForce_[i] << " ";
         Info << endl;
    }
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
        
        // Debug
        /*
        if (Pstream::master() && i==0)
        {
             Info<< "Mode 0: Freq=" << freq << " F_this=" << F_this << " F_last=" << F_last << endl;
        }
        */

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

        /*
        if (Pstream::master() && i==0)
        {
             Info << "Mode 0 Step:"
                  << " dt=" << dt
                  << " theta=" << theta_
                  << " F_this=" << F_this 
                  << " load=" << load 
                  << " effK=" << effectiveK 
                  << " disTheta=" << disTheta
                  << " dis=" << dis 
                  << " acc=" << acc << endl;
        }
        */
        
        modeState_[i] = vector(dis, vel, acc);
    }
}

void fastDynamicFvMesh::writeDiagnostics() const
{
    if (!Pstream::master())
    {
        return;
    }

    std::ofstream diagFile;
    bool writeHeader = false;

    std::ifstream check("modal_diagnostics.csv");
    if (!check.good())
    {
        writeHeader = true;
    }
    check.close();

    diagFile.open("modal_diagnostics.csv", std::ios::app);

    if (!diagFile.good())
    {
        WarningInFunction << "Unable to open modal_diagnostics.csv for writing." << endl;
        return;
    }

    if (writeHeader)
    {
        diagFile << "Time";
        for (label i=0; i<nMode_; ++i)
        {
            diagFile << ",Force_" << (i+1);
        }

        for (label i=0; i<nMode_; ++i)
        {
            diagFile << ",Disp_" << (i+1)
                     << ",Vel_" << (i+1)
                     << ",Acc_" << (i+1);
        }

        diagFile << "\n";
    }

    diagFile << std::setprecision(12);
    diagFile << this->time().value();

    for (label i=0; i<nMode_; ++i)
    {
        diagFile << "," << modeForce_[i];
    }

    for (label i=0; i<nMode_; ++i)
    {
        diagFile << "," << modeState_[i].x()
                 << "," << modeState_[i].y()
                 << "," << modeState_[i].z();
    }

    diagFile << "\n";
}

bool fastDynamicFvMesh::update()
{
    if (Pstream::master()) Info<< "DEBUG: Starting update()" << endl;

    // 1. Calculate time step
    scalar dt = this->time().deltaTValue();
    
    // Store previous state
    modeForce0_ = modeForce_;
    modeState0_ = modeState_;

    // 2. Calculate forces
    if (Pstream::master()) Info<< "DEBUG: Calling calcModalForces()" << endl;
    calcModalForces();
    if (Pstream::master()) Info<< "DEBUG: Finished calcModalForces()" << endl;

    if (startupStepCount_ < 2)
    {
        if (Pstream::master())
        {
            Info<< "Initializing Fast Dynamic Mesh state (startup step "
                << (startupStepCount_ + 1) << " of 2)..." << endl;
        }

        for (label i=0; i<nMode_; ++i)
        {
            scalar F = modeForce_[i];
            scalar v = 0.0;
            if (i < initVelocity_.size()) v = initVelocity_[i];

            // Legacy Fluent UDF initialisation:
            // D0 = 0, V0 = initVelocity, A0 = F0
            scalar a = F;

            modeState_[i] = vector(0.0, v, a);
            modeState0_[i] = modeState_[i];
            modeForce0_[i] = modeForce_[i];
        }

        ++startupStepCount_;
        writeDiagnostics();

        if (Pstream::master())
        {
            Info<< "DEBUG: Skipping mesh motion during startup initialisation" << endl;
        }

        return true;
    }

    // 3. Solve dynamics
    if (nMode_ > 0)
    {
        if (Pstream::master()) Info<< "DEBUG: Calling solveStructuralDynamics()" << endl;
        solveStructuralDynamics(dt);
        if (Pstream::master()) Info<< "DEBUG: Finished solveStructuralDynamics()" << endl;
    }

    // 4. Update mesh points
    if (Pstream::master()) Info<< "DEBUG: Updating mesh points" << endl;
    pointField newPoints = this->points(); 
 

    // Debug output:
    if (Pstream::master() && nMode_ > 0)
    {
        Info<< "Time: " << this->time().value() << " First Mode Disp: " << modeState_[0].x() << endl;
    }

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
    writeDiagnostics();

    return true;
}

} // End namespace Foam
