/*----------------------------------------------------------------***********\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Fast Dynamic Mesh method based on modal superposition.

\*---------------------------------------------------------------------------*/

#ifndef fastDynamicFvMesh_H
#define fastDynamicFvMesh_H

#include "dynamicFvMesh.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class fastDynamicFvMesh Declaration
\*---------------------------------------------------------------------------*/

class fastDynamicFvMesh
:
    public dynamicFvMesh
{
    // Private Data

        //- Number of modes
        label nMode_;

        //- Modal frequencies
        scalarField modeFreq_;

        //- Modal damping ratios
        scalarField modeDamp_;

        //- Modal stiffness
        scalarField modeStiff_;

        //- Modal force (current time step)
        scalarField modeForce_;

        //- Pressure contribution to modal force
        scalarField modePressureForce_;

        //- Viscous wall-shear contribution to modal force
        scalarField modeShearForce_;

        //- Modal force (previous time step)
        scalarField modeForce0_;

        //- Modal displacement (current time step)
        //  Layout: [disp, vel, acc] per mode
        List<vector> modeState_;

        //- Modal displacement (previous time step)
        List<vector> modeState0_;

        //- Initial modal velocity
        scalarField initVelocity_;

        //- Relaxed modal displacement actually imposed on the mesh
        scalarField appliedModeDisp_;

        //- Mode shapes (displacement field per mode)
        //  Stores full field for simplicity
        List<vectorField> modeShapes_;

        //- Wilson-Theta parameter
        scalar theta_;

        //- Point-mapping tolerance used for legacy CSV node matching
        scalar mappingTolerance_;

        //- Relaxation applied to structural modal displacement before moving mesh
        scalar couplingRelaxation_;

        //- Pressure field used for modal loading
        word pressureFieldName_;

        //- Reference density for kinematic pressure fields
        scalar rhoRef_;

        //- FSI patch names
        wordRes fsiPatches_;

        //- Write per-face modal-force diagnostics on the FSI interface
        bool writeFaceDiagnostics_;

        //- Zero-based mode index for face diagnostics
        label faceDiagnosticsMode_;

        //- Number of startup calls already treated as Wilson-theta initialisation
        label startupStepCount_;

        //- Last physical time index already processed by update()
        label lastUpdateTimeIndex_;


    // Private Member Functions

        //- Read modal properties from dictionary and files
        void readControls();

        //- Read mode shapes from external files
        void readModeShapes();

        //- Read legacy runtime parameters from FluidPara.csv
        void readLegacyParameters(const fileName& modeDir);

        //- Calculate modal forces from fluid pressure and shear stress
        void calcModalForces();

        //- Return patch density used for shear calculations
        tmp<scalarField> patchDensity(const label patchi, scalar defaultRho) const;

        //- Deviatoric effective stress tensor for wall shear evaluation
        tmp<symmTensorField> devRhoReff
        (
            const tensorField& gradUp,
            const label patchi,
            scalar defaultRho
        ) const;

        //- Solve structural dynamics using Wilson-Theta method
        void solveStructuralDynamics(scalar dt);

        //- Write modal force/displacement diagnostics
        void writeDiagnostics() const;

        //- Disallow default bitwise copy construct
        fastDynamicFvMesh(const fastDynamicFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const fastDynamicFvMesh&);


public:

    //- Runtime type information
    TypeName("fastDynamicFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit fastDynamicFvMesh(const IOobject& io);


    //- Destructor
    virtual ~fastDynamicFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
