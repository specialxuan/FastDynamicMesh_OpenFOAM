# Fast Dynamic Mesh for OpenFOAM

This library implements a fast dynamic mesh method based on modal superposition, ported from an ANSYS Fluent UDF.

## Usage

1. Compile the library:
   wmake libso

2. Add the library to your `controlDict`:
   libs ( "libfastDynamicFvMesh.so" );

3. Configure `constant/dynamicMeshDict`:
   dynamicFvMesh fastDynamicFvMesh;

   fastDynamicFvMeshCoeffs
   {
    theta 1.4; // Wilson-Theta parameter
    fsiPatches ( "wall" ); // List of patches where fluid forces are calculated
    mappingTolerance 4e-6; // Optional CSV-to-mesh node matching tolerance
   }

4. Fluid loading
   Modal forces on `fsiPatches` now include both pressure and viscous wall-shear contributions, projected onto the averaged face mode shape to mirror the legacy Fluent force assembly as closely as practical. Ensure the velocity field `U` is available so wall shear can be evaluated.

5. Input Files
   Place the mode shape files in a `mode/` directory in your case root:
   - `mode/FluidNodeCoor.csv`
   - `mode/FluidNodeDisp1.csv`, `mode/FluidNodeDisp2.csv`, etc.
   - Optional: `mode/FluidPara.csv` to provide legacy initial modal velocities.

   `FluidPara.csv` legacy zone IDs are read for compatibility and logged, but the
   OpenFOAM implementation still selects FSI surfaces by `fsiPatches`.
