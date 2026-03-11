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
     couplingRelaxation 0.1; // Optional relaxation on displacement applied to the fluid mesh
     pressureField p; // Optional pressure field name, defaults to p
     rhoRef 1000; // Required for kinematic-pressure fields unless rho is available
    }

4. Fluid loading
   Modal forces on `fsiPatches` now include both pressure and viscous wall-shear contributions, projected onto the averaged face mode shape to mirror the legacy Fluent force assembly as closely as practical. Ensure the velocity field `U` is available so wall shear can be evaluated.
   Pressure loading is read from `pressureField` (default `p`). If that field has pressure dimensions, it is used directly. If it has kinematic-pressure dimensions, the mesh requires either `rhoRef`, a `rho` entry in `transportProperties`, or a `rho` volScalarField so the pressure can be converted to force consistently.
   `couplingRelaxation` relaxes the structural displacement before it is imposed on the fluid mesh. Values below `1` add coupling-side damping to help stabilize loosely coupled FSI runs without changing the structural solve itself.
   The run-time `modal_diagnostics.csv` reports total modal force together with separate `PressureForce_i` and `ShearForce_i` columns to help debug force decomposition. It also records the structural modal displacement/state and the relaxed `AppliedDisp_i` actually sent to the mesh.
   For face-level debugging on one mode, enable `writeFaceDiagnostics yes;` and set `faceDiagnosticsMode` (1-based) in `fastDynamicFvMeshCoeffs`. Parallel runs write per-processor CSV files named like `faceDiagnostics_mode4_proc0.csv`.

5. Input Files
   Place the mode shape files in a `mode/` directory in your case root:
   - `mode/FluidNodeCoor.csv`
   - `mode/FluidNodeDisp1.csv`, `mode/FluidNodeDisp2.csv`, etc.
   - Optional: `mode/FluidPara.csv` to provide legacy initial modal velocities.

   `FluidPara.csv` legacy zone IDs are read for compatibility and logged, but the
   OpenFOAM implementation still selects FSI surfaces by `fsiPatches`.
