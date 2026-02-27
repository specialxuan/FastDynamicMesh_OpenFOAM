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

        //- Modal force (previous time step)
        scalarField modeForce0_;

        //- Modal displacement (current time step)
        //  Layout: [disp, vel, acc] per mode
        List<vector> modeState_;

        //- Modal displacement (previous time step)
        List<vector> modeState0_;

        //- Initial modal velocity
        scalarField initVelocity_;

        //- Mode shapes (displacement field per mode)
        //  Stores full field for simplicity
        List<vectorField> modeShapes_;

        //- Wilson-Theta parameter
        scalar theta_;

        //- FSI patch names
        wordRes fsiPatches_;


    // Private Member Functions

        //- Read modal properties from dictionary and files
        void readControls();

        //- Read mode shapes from external files
        void readModeShapes();

        //- Calculate modal forces from fluid pressure and shear stress
        void calcModalForces();

        //- Solve structural dynamics using Wilson-Theta method
        void solveStructuralDynamics(scalar dt);

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
