Fast Dynamic Mesh for OpenFOAM
==============================

This library implements a fast dynamic mesh method based on modal superposition, ported from an ANSYS Fluent UDF.

Usage
-----
1. Compile the library:
   wmake libso

2. Add the library to your `controlDict`:
   libs ( "libfastDynamicFvMesh.so" );

3. Configure `constant/dynamicMeshDict`:
   dynamicFvMesh   fastDynamicFvMesh;
   
   fastDynamicFvMeshCoeffs
   {
       theta       1.4;            // Wilson-Theta parameter
       fsiPatches  ( "wall" );     // List of patches where fluid forces are calculated
   }

4. Input Files
   Place the mode shape files in a `mode/` directory in your case root:
   - `mode/FluidNodeCoor.csv`
   - `mode/FluidNodeDisp1.csv`, `mode/FluidNodeDisp2.csv`, etc.
