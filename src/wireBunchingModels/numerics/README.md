# beamFoam numerics

Shared mathematical and utility code for the beam library. Check here before
writing a new helper: extend an existing function rather than adding a second
version of it. When you add, rename or remove a helper, update this list in the
same change.

## Rotations

- `pseudoVector/`
  - `rotationMatrix(...)`: rotation matrix from a rotation vector (Rodrigues'
    formula). Overloads for `vector`, `vectorField`, `volVectorField` and
    `surfaceVectorField`.
  - `pseudoVector(...)`: rotation vector from a rotation matrix (Spurrier's
    algorithm via quaternions). Safe at a half-turn; the angle is in
    `[0, pi]`.
  - `rotationTangentOperator(...)`: tangent operator (left Jacobian, dexp) of
    a rotation vector. `T` maps a change in the rotation vector to the spatial
    rotation and `T.T()` to the material rotation.
  - `interpolateRotationMatrix(...)`: interpolation of rotation matrices
    between cell and face centres.
  - The header lists which functions the solver uses and which are kept for
    reference only.
- `spinTensor/`
  - `spinTensor(...)`: skew-symmetric tensor from an axial vector.
  - `axialVector(...)`: axial vector from a skew-symmetric tensor.
  - Overloads for single values, fields, `vol` and `surface` fields.

## Geometry

- `HermiteSpline/`: cubic Hermite spline through points and tangents. Gives
  positions, derivatives, curvatures and lengths along the beam centreline.
- `beamHelperFunctions/`
  - `sharedFace(c0, c1)`: the face shared by two cells.
  - `neiPatchIndex(c, mesh)`: the patch of a cell's boundary face.

## Linear algebra

- `BlockEigenSolvers/BlockEigenSolverOF/`: solves the 6x6 block-coupled
  Newton system of the beam with Eigen.
- `denseMatrixHelperFunctions/`: `magSqr(...)` of a field of
  `scalarRectangularMatrix`.
- `FieldSumOp/`: element-wise sum of two fields, for use as the operator in
  `reduce` and similar calls.

## Input types

- `Tuples/`: `Tuple2` types (`labelScalar`, `labelScalarVector`, ...) read
  from beam dictionaries.
- `Tuples/timeVaryingVector/`: a vector that is constant or interpolated from
  a time series.

## Fluid coupling and parallel checks

- `samplingFluid/`
  - `getFluidVelocity(...)`: samples the fluid velocity around the beam for
    beam-fluid coupling (used by `almDrag`).
  - `upstreamSampling.H`: upstream sampling points for the actuator-line
    model.
- `beamParallelChecks/`: `checkBeamOnMasterProcessor(...)` stops with an
  error if a beam region that must stay on the master processor has been
  decomposed.

## Elsewhere in the library

- `beamModels/beamModel/lookupBeamModel.H`: find the beam model of a mesh
  region.
- `beamModels/beamModel/beamSolverCounters.H`: per-time-step counts of
  `evolve()` calls and Newton iterations.
- `beamMomentumContribution/beamMomentumContribution.H`: `checkSingleBeam`,
  `checkCircularSection` and `unitTangents` for contributions.
- OpenFOAM's own `src` (for example `meshTools`, `finiteVolume`,
  `OpenFOAM/primitives`) often has what you need already.
