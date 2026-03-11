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
    couplingRelaxation_(1.0),
    pressureFieldName_("p"),
    rhoRef_(-1.0),
    writeFaceDiagnostics_(false),
    faceDiagnosticsMode_(-1),
    startupStepCount_(0),
    lastUpdateTimeIndex_(-1)
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
    
    if (!dynamicMeshDict.found(typeName + "Coeffs"))
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Missing required sub-dictionary '" << typeName << "Coeffs' in "
            << dynamicMeshDict.objectPath() << nl
            << "Add constant/dynamicMeshDict entries for " << typeName << '.'
            << exit(FatalIOError);
    }

    const dictionary& fdmDict = dynamicMeshDict.subDict(typeName + "Coeffs");

    fdmDict.readIfPresent("theta", theta_);

    if (!fdmDict.found("fsiPatches"))
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Missing required entry 'fsiPatches' in sub-dictionary '"
            << typeName << "Coeffs' of " << dynamicMeshDict.objectPath()
            << exit(FatalIOError);
    }

    fdmDict.lookup("fsiPatches") >> fsiPatches_;

    if (!fsiPatches_.size())
    {
        FatalErrorInFunction
            << "At least one FSI patch name must be provided in "
            << dynamicMeshDict.objectPath() << exit(FatalError);
    }

    fdmDict.readIfPresent("mappingTolerance", mappingTolerance_);
    fdmDict.readIfPresent("couplingRelaxation", couplingRelaxation_);
    fdmDict.readIfPresent("pressureField", pressureFieldName_);
    const bool foundRhoRef = fdmDict.readIfPresent("rhoRef", rhoRef_);

    if (couplingRelaxation_ <= 0 || couplingRelaxation_ > 1)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'couplingRelaxation' must be in the range (0, 1] in "
            << "sub-dictionary '" << typeName << "Coeffs' of "
            << dynamicMeshDict.objectPath() << exit(FatalIOError);
    }

    if (foundRhoRef && rhoRef_ <= 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "Entry 'rhoRef' must be positive in sub-dictionary '"
            << typeName << "Coeffs' of " << dynamicMeshDict.objectPath()
            << exit(FatalIOError);
    }

    label faceDiagnosticsMode = -1;
    const bool foundFaceDiagnosticsMode =
        fdmDict.readIfPresent("faceDiagnosticsMode", faceDiagnosticsMode);

    fdmDict.readIfPresent("writeFaceDiagnostics", writeFaceDiagnostics_);

    if (foundFaceDiagnosticsMode)
    {
        writeFaceDiagnostics_ = true;
        faceDiagnosticsMode_ = faceDiagnosticsMode - 1;
    }

    if (writeFaceDiagnostics_ && faceDiagnosticsMode_ < 0)
    {
        FatalIOErrorInFunction(dynamicMeshDict)
            << "When 'writeFaceDiagnostics' is enabled, provide a positive "
            << "'faceDiagnosticsMode' entry in sub-dictionary '"
            << typeName << "Coeffs' of " << dynamicMeshDict.objectPath()
            << exit(FatalIOError);
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
            FatalErrorInFunction
                << "Cannot open required mode coordinate file " << coorPath << nl
                << "Provide mode/FluidNodeCoor.csv in the case root before "
                << "running fastDynamicFvMesh." << exit(FatalError);
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
            if (!(ss >> dummy >> nNode >> nMode) || nNode <= 0 || nMode <= 0)
            {
                FatalErrorInFunction
                    << "Invalid FluidNodeCoor.csv header in " << coorPath << nl
                    << "Expected three comma-separated values with positive "
                    << "node and mode counts." << exit(FatalError);
            }
            
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
                if (!std::getline(file, line))
                {
                    FatalErrorInFunction
                        << "Unexpected end of file while reading node " << i
                        << " from " << coorPath << exit(FatalError);
                }
                std::replace(line.begin(), line.end(), ',', ' ');
                std::stringstream ss2(line);
                if (!(ss2 >> x >> y >> z))
                {
                    FatalErrorInFunction
                        << "Invalid coordinate entry for node " << i << " in "
                        << coorPath << exit(FatalError);
                }
                
                csvPoints[i] = point(x, y, z);
            }
            
            // Read mode shapes
            for (label m=0; m<nMode_; ++m)
            {
                fileName shapePath = modeDir/("FluidNodeDisp" + std::to_string(m+1) + ".csv");
                std::ifstream mFile(shapePath);

                if (!mFile.good())
                {
                    FatalErrorInFunction
                        << "Cannot open required mode shape file " << shapePath
                        << exit(FatalError);
                }

                string line;
                if (!std::getline(mFile, line))
                {
                    FatalErrorInFunction
                        << "Missing frequency line in " << shapePath
                        << exit(FatalError);
                }

                std::replace(line.begin(), line.end(), ',', ' ');
                std::stringstream ss(line);
                scalar freq = 0.0;
                scalar fileNodeCount = 0.0;
                scalar fileModeCount = 0.0;

                if (!(ss >> freq))
                {
                    FatalErrorInFunction
                        << "Failed to read frequency from " << shapePath
                        << exit(FatalError);
                }

                if ((ss >> fileNodeCount) && (ss >> fileModeCount))
                {
                    if
                    (
                        label(fileNodeCount) != nCsvNodes
                     || label(fileModeCount) != nMode_
                    )
                    {
                        FatalErrorInFunction
                            << "Mode file " << shapePath
                            << " reports " << label(fileNodeCount)
                            << " nodes and " << label(fileModeCount)
                            << " modes, but FluidNodeCoor.csv reports "
                            << nCsvNodes << " nodes and " << nMode_
                            << " modes." << exit(FatalError);
                    }
                }

                modeFreq_[m] = freq;
                Info<< "  Mode " << m << " Freq: " << freq << endl;

                label dataRow = 0;

                while (dataRow < nCsvNodes)
                {
                    if (!std::getline(mFile, line))
                    {
                        FatalErrorInFunction
                            << "Unexpected end of file while reading node "
                            << dataRow << " from " << shapePath
                            << exit(FatalError);
                    }

                    std::replace(line.begin(), line.end(), ',', ' ');
                    std::stringstream ss2(line);

                    scalar dx = 0.0;
                    scalar dy = 0.0;
                    scalar dz = 0.0;

                    if (!(ss2 >> dx >> dy >> dz))
                    {
                        if (dataRow == 0)
                        {
                            continue;
                        }

                        FatalErrorInFunction
                            << "Invalid displacement entry for node "
                            << dataRow << " in " << shapePath
                            << exit(FatalError);
                    }

                    csvShapes[m][dataRow] = vector(dx, dy, dz);
                    ++dataRow;
                }
            }
        }
    }

    // Broadcast sizes
    Pstream::broadcast(nMode_);
    Pstream::broadcast(nCsvNodes);
    
    // Resize local arrays
    if (!Pstream::master()) modeFreq_.setSize(nMode_); // Only slaves need resize now

    if (nMode_ <= 0)
    {
        FatalErrorInFunction
            << "No fluid modes were loaded from " << modeDir << nl
            << "Ensure mode/FluidNodeCoor.csv and mode/FluidNodeDisp*.csv are "
            << "present and valid." << exit(FatalError);
    }

    if (writeFaceDiagnostics_ && faceDiagnosticsMode_ >= nMode_)
    {
        FatalErrorInFunction
            << "Requested faceDiagnosticsMode " << (faceDiagnosticsMode_ + 1)
            << " but only " << nMode_ << " modes were loaded."
            << exit(FatalError);
    }

    modeForce_.setSize(nMode_, 0.0);
    modePressureForce_.setSize(nMode_, 0.0);
    modeShearForce_.setSize(nMode_, 0.0);
    modeForce0_.setSize(nMode_, 0.0);
    modeState_.setSize(nMode_, vector::zero);
    modeState0_.setSize(nMode_, vector::zero);
    initVelocity_.setSize(nMode_, 0.0);
    appliedModeDisp_.setSize(nMode_, 0.0);

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

    const label totalMapped = returnReduce(mappedCount, sumOp<label>());
    const label totalPoints = returnReduce(label(localPoints.size()), sumOp<label>());

    if (Pstream::master())
    {
        Info<< "Mapped " << totalMapped << " points out of "
            << totalPoints << " mesh points." << endl;
    }

    if (totalMapped == 0)
    {
        FatalErrorInFunction
            << "Mapped 0 mesh points from " << modeDir << nl
            << "Check that the mode CSV coordinates correspond to the current "
            << "mesh and that mappingTolerance is large enough." << exit(FatalError);
    }

    if (Pstream::master() && totalMapped != totalPoints)
    {
        WarningInFunction
            << "Mapped " << totalMapped << " of " << totalPoints
            << " mesh points. Unmapped points will remain stationary."
            << endl;
    }
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
    modePressureForce_ = 0.0;
    modeShearForce_ = 0.0;

    if (!this->foundObject<volScalarField>(pressureFieldName_))
    {
        FatalErrorInFunction
            << "Required pressure field '" << pressureFieldName_
            << "' was not found. Configure 'pressureField' in dynamicMeshDict "
            << "or supply the field before fastDynamicFvMesh::update()."
            << exit(FatalError);
    }

    const volScalarField& p = this->lookupObject<volScalarField>(pressureFieldName_);
    const dimensionSet kinematicPressureDims(dimPressure/dimDensity);
    const bool dimensionalPressure = (p.dimensions() == dimPressure);
    const bool kinematicPressure = (p.dimensions() == kinematicPressureDims);

    if (!dimensionalPressure && !kinematicPressure)
    {
        FatalErrorInFunction
            << "Pressure field '" << pressureFieldName_ << "' has dimensions "
            << p.dimensions() << ". Expected either " << dimPressure
            << " (pressure) or " << kinematicPressureDims
            << " (kinematic pressure)." << exit(FatalError);
    }

    scalar defaultRho = rhoRef_;

    if (kinematicPressure && defaultRho <= 0)
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

        if (transportProperties.found("rho"))
        {
            transportProperties.lookup("rho") >> defaultRho;
        }
    }

    if (kinematicPressure && defaultRho <= 0 && !this->foundObject<volScalarField>("rho"))
    {
        FatalErrorInFunction
            << "Pressure field '" << pressureFieldName_
            << "' is kinematic, but no density information was found." << nl
            << "Provide 'rhoRef' in dynamicMeshDict, a 'rho' entry in "
            << "constant/transportProperties, or a volScalarField named 'rho'."
            << exit(FatalError);
    }

    std::ofstream faceDiagFile;
    bool writeFaceDiagHeader = false;
    bool haveFaceDiagFile = false;

    if (writeFaceDiagnostics_)
    {
        fileName diagnosticsRoot = this->time().path();

        if (Pstream::parRun())
        {
            diagnosticsRoot = this->time().path().path();
        }

        const word modeLabel("mode" + Foam::name(faceDiagnosticsMode_ + 1));
        fileName faceDiagnosticsPath =
            diagnosticsRoot/
            ("faceDiagnostics_" + modeLabel + "_proc"
           + Foam::name(Pstream::myProcNo()) + ".csv");

        std::ifstream check(faceDiagnosticsPath.c_str());
        writeFaceDiagHeader = !check.good();
        check.close();

        faceDiagFile.open(faceDiagnosticsPath.c_str(), std::ios::app);

        if (!faceDiagFile.good())
        {
            WarningInFunction
                << "Unable to open " << faceDiagnosticsPath
                << " for face-level diagnostics." << endl;
        }
        else
        {
            haveFaceDiagFile = true;
            faceDiagFile << std::setprecision(12);

            if (writeFaceDiagHeader)
            {
                faceDiagFile
                    << "Processor,Time,Patch,PatchFace,MeshFace,"
                    << "Cx,Cy,Cz,"
                    << "AreaX,AreaY,AreaZ,"
                    << "ShapeX,ShapeY,ShapeZ,"
                    << "PressureForceX,PressureForceY,PressureForceZ,"
                    << "ShearForceX,ShearForceY,ShearForceZ,"
                    << "PressureContribution,ShearContribution,TotalContribution\n";
            }
        }
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

    // Debug p range
    scalar pMin = gMin(p);
    scalar pMax = gMax(p);
    if (Pstream::master())
    {
        Info << "DEBUG: " << pressureFieldName_ << " range [" << pMin << ", "
             << pMax << "]" << endl;

        if (kinematicPressure)
        {
            Info << "DEBUG: pressure field is kinematic; density comes from "
                 << (this->foundObject<volScalarField>("rho") ? "field rho" : "rhoRef/transportProperties")
                 << endl;
        }
    }

    // Iterate over FSI patches
    forAll(fsiPatches_, i)
    {
        label patchID = this->boundaryMesh().findPatchID(fsiPatches_[i]);
        if (patchID == -1)
        {
            FatalErrorInFunction
                << "Configured FSI patch '" << fsiPatches_[i]
                << "' was not found in boundaryMesh." << exit(FatalError);
        }

        const polyPatch& pp = this->boundaryMesh()[patchID];
        const fvPatchScalarField& pPatch = p.boundaryField()[patchID];
        const vectorField faceAreas(pp.faceAreas());
        const vectorField faceCentres(pp.faceCentres());

        tmp<scalarField> tPressureScale
        (
            new scalarField(pp.size(), 1.0)
        );

        if (kinematicPressure)
        {
            tPressureScale = patchDensity(patchID, defaultRho);
        }

        const scalarField& pressureScale = tPressureScale();

        const symmTensorField* devStressPtr = nullptr;
        tmp<symmTensorField> tDevStress;

        if (gradUPtr)
        {
            const tensorField& gradPatch = gradUPtr->boundaryField()[patchID];
            tDevStress = devRhoReff(gradPatch, patchID, defaultRho > 0 ? defaultRho : 1.0);
            devStressPtr = &tDevStress();
        }

        // Loop over faces
        forAll(pp, faceI)
        {
            const vector& areaVec = faceAreas[faceI];
            const vector pressureForce = pressureScale[faceI] * pPatch[faceI] * areaVec;
            vector shearForce = vector::zero;

            if (devStressPtr)
            {
                shearForce = areaVec & (*devStressPtr)[faceI];
                const scalar areaMag = mag(areaVec);

                if (areaMag > VSMALL)
                {
                    const vector faceNormal = areaVec/areaMag;
                    shearForce -= faceNormal*(faceNormal & shearForce);
                }
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

                const scalar pressureModeForce = pressureForce & shapeFace;
                const scalar shearModeForce = shearForce & shapeFace;

                modePressureForce_[m] += pressureModeForce;
                modeShearForce_[m] += shearModeForce;
                modeForce_[m] += pressureModeForce + shearModeForce;

                if (haveFaceDiagFile && m == faceDiagnosticsMode_)
                {
                    const point& faceCentre = faceCentres[faceI];

                    faceDiagFile
                        << Pstream::myProcNo() << ','
                        << this->time().value() << ','
                        << pp.name() << ','
                        << faceI << ','
                        << (pp.start() + faceI) << ','
                        << faceCentre.x() << ','
                        << faceCentre.y() << ','
                        << faceCentre.z() << ','
                        << areaVec.x() << ','
                        << areaVec.y() << ','
                        << areaVec.z() << ','
                        << shapeFace.x() << ','
                        << shapeFace.y() << ','
                        << shapeFace.z() << ','
                        << pressureForce.x() << ','
                        << pressureForce.y() << ','
                        << pressureForce.z() << ','
                        << shearForce.x() << ','
                        << shearForce.y() << ','
                        << shearForce.z() << ','
                        << pressureModeForce << ','
                        << shearModeForce << ','
                        << (pressureModeForce + shearModeForce) << '\n';
                }
            }
        }
    }

    // Parallel reduction
    Pstream::listCombineGather(modeForce_, plusEqOp<scalar>());
    Pstream::listCombineGather(modePressureForce_, plusEqOp<scalar>());
    Pstream::listCombineGather(modeShearForce_, plusEqOp<scalar>());
    Pstream::broadcast(modeForce_);
    Pstream::broadcast(modePressureForce_);
    Pstream::broadcast(modeShearForce_);

    if (Pstream::master())
    {
        Info << "DEBUG: Mode Forces p/shear/total (0-4): ";
        for (label i=0; i<min(5, nMode_); ++i)
        {
            Info << '['
                 << modePressureForce_[i] << ','
                 << modeShearForce_[i] << ','
                 << modeForce_[i] << "] ";
        }
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

    fileName diagnosticsPath = this->time().path()/"modal_diagnostics.csv";

    if (Pstream::parRun())
    {
        diagnosticsPath = this->time().path().path()/"modal_diagnostics.csv";
    }

    std::ofstream diagFile;
    bool writeHeader = false;

    std::ifstream check(diagnosticsPath.c_str());
    if (!check.good())
    {
        writeHeader = true;
    }
    check.close();

    diagFile.open(diagnosticsPath.c_str(), std::ios::app);

    if (!diagFile.good())
    {
        WarningInFunction << "Unable to open " << diagnosticsPath
            << " for writing." << endl;
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
            diagFile << ",PressureForce_" << (i+1);
        }

        for (label i=0; i<nMode_; ++i)
        {
            diagFile << ",ShearForce_" << (i+1);
        }

        for (label i=0; i<nMode_; ++i)
        {
            diagFile << ",Disp_" << (i+1)
                     << ",Vel_" << (i+1)
                     << ",Acc_" << (i+1)
                     << ",AppliedDisp_" << (i+1);
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
        diagFile << "," << modePressureForce_[i];
    }

    for (label i=0; i<nMode_; ++i)
    {
        diagFile << "," << modeShearForce_[i];
    }

    for (label i=0; i<nMode_; ++i)
    {
        diagFile << "," << modeState_[i].x()
                 << "," << modeState_[i].y()
                 << "," << modeState_[i].z()
                 << "," << appliedModeDisp_[i];
    }

    diagFile << "\n";
}

bool fastDynamicFvMesh::update()
{
    if (Pstream::master()) Info<< "DEBUG: Starting update()" << endl;

    const label currentTimeIndex = this->time().timeIndex();
    if (currentTimeIndex == lastUpdateTimeIndex_)
    {
        return false;
    }
    lastUpdateTimeIndex_ = currentTimeIndex;

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
            appliedModeDisp_[i] = 0.0;
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
        Info<< "Time: " << this->time().value()
            << " First Mode Disp: " << modeState_[0].x()
            << " Applied: " << appliedModeDisp_[0] << endl;
    }

    for (label m=0; m<nMode_; ++m)
    {
        const vectorField& shape = modeShapes_[m];

        const scalar appliedDisp0 = appliedModeDisp_[m];
        const scalar appliedDisp =
            appliedDisp0
          + couplingRelaxation_*(modeState_[m].x() - appliedDisp0);

        // Incremental displacement actually imposed on the mesh
        scalar dDisp = appliedDisp - appliedDisp0;
        appliedModeDisp_[m] = appliedDisp;

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
