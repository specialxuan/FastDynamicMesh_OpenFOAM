# Fast Dynamic Mesh Adapter - Session Summary

## Overview

This session focused on fixing the integration of the Fast Dynamic Mesh (FDM) method into OpenFOAM (v2412) and validating it with the `validationCase`.

## Changes Made

### 1. Library Fixes (`src/FastDynamic/fastDynamicFvMesh/`)

- **Fixed `readControls()` in `fastDynamicFvMesh.C`**:
  - Replaced `lookupObject<dictionary>("dynamicMeshDict")` which failed because the dictionary wasn't registered.
  - Implemented direct `IOdictionary` reading of `constant/dynamicMeshDict` to ensure parameters like `theta` and `fsiPatches` are loaded correctly.
- **Fixed Parallel Broadcast Bug**:
  - **Issue**: `csvShapes` array was not resized on slave processors before `Pstream::broadcast`, causing a segmentation fault in parallel runs.
  - **Fix**: Added `csvShapes.setSize(nMode_)` before the broadcast loop in `readModeShapes()`.
- **Fixed Includes**: Added missing includes (`Pstream.H`) to `fastDynamicFvMesh.C`.
- **Fixed Path Resolution**:
  - Changed file reading to use absolute paths (`this->time().path()/"mode"/...`) and added fallback to parent directory (`path().path()/"mode"`) for parallel execution.
- **Improved CSV Parsing**:
  - Rewrote `readModeShapes` to robustly handle trailing commas and whitespace in CSV files, which previously caused read failures ("Mapped 0 points").
- **Relaxed Point Mapping Tolerance**:
  - Changed distance tolerance `tolSqr` from `1e-12` to `1e-8` to handle floating point discrepancies between CSV coordinates and OpenFOAM mesh nodes. This resolved the critical "Mapped 0 points" failure.

### 2. Validation Case Fixes (`run/validationCase/`)

- **`0/U` (Velocity Boundary Condition)**:
  - Implemented parabolic inlet using `timeVaryingMappedFixedValue` to bypass `codedFixedValue` security restrictions (root user compilation error).
  - Created `generate_inlet.py` to pre-calculate the boundary data in `constant/boundaryData`.
- **`system/fvSolution`**:
  - Added `pcorr` (pressure correction) and `pcorrFinal` solvers.
  - **Reason**: `pimpleFoam` with dynamic mesh motion requires a pressure correction step to satisfy continuity on the moving mesh, which was missing from the solver settings.
- **`run_parallel.sh`**:
  - Created a robust automation script to handle cleanup, inlet generation, decomposition, parallel solving, and reconstruction (including mesh reconstruction).

### 3. Diagnostics and Scaling

- **Fixed Force Scaling**:
  - **Issue**: Forces were negligible (~1e-3 N) in `pimpleFoam` (incompressible solver) because it uses kinematic pressure ($p/\rho$).
  - **Fix**: Added logic to read density `rho` from `constant/transportProperties` and scale forces by density.
  - **Result**: Modal forces and displacements now match expected physical magnitudes.
- **Diagnostics Output**:
  - Implemented `modal_diagnostics.csv` output to track forces and displacements at each time step.

## Verification

- **Compilation**: The `fastDynamicFvMesh` library was successfully compiled using `wmake`.
- **Execution**: `pimpleFoam` runs successfully in parallel.
- **Mesh Motion**: Verified that mesh points move (checked via `diff` of `polyMesh/points` and diagnostics file).
- **Mapping**: Confirmed successful mapping of ~25k points (was 0).

## How to Run

1. Source OpenFOAM environment (if not already done).
2. Compile the library (if needed):

   ```bash
   cd src/FastDynamic/fastDynamicFvMesh
   wmake
   ```

3. Run the validation case:

   ```bash
   cd ~/Workspace/validationCase
   ./run_parallel.sh
   ```

### 4. Stability Fixes (Run-Time)

- **Crash at t=0.441s**:
  - **Issue**: Divergence (Courant number > 22) leading to Floating Point Exception.
  - **Fix**:
    - **Adaptive Time Stepping**: Enabled in `controlDict` (`adjustTimeStep yes`, `maxCo 1.0`, `maxDeltaT 0.005`).
    - **Solver Robustness**: Increased `nOuterCorrectors` (3) and `nCorrectors` (3) in `fvSolution`.
    - **Schemes**: Switched `div(phi,U)` from `Gauss LUST` to `Gauss upwind` temporarily to stabilize the run.
