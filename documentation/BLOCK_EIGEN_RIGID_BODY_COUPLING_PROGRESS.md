# BlockEigen Rigid-Body Coupling Progress

Date: 2026-05-07

This note summarises the work done to move the `rigidBodyAndBeam` case from
using the original sixDoF beam-restraint force path toward using the rigid-body
motion solved inside `BlockEigenSolverOF`.

## Objective

The original coupled case let `sixDoFRigidBodyMotionFvBeam` update the rigid-body
state, while beamFoam supplied a restraint force back to that solver. The newer
work added six rigid-body rows to the beamFoam `BlockEigenSolverOF`, so the next
goal was to make that BlockEigen rigid-body solution usable by the moorFV
rigid-body solver.

The immediate milestone was:

1. Have `BlockEigenSolverOF` solve and expose the full rigid-body state update.
2. Pass that solution back out through the beam model interface.
3. Let the moorFV `finiteVolumeBeam` restraint optionally make the BlockEigen
   solution authoritative for the rigid body.
4. Avoid applying the beam restraint force twice.
5. Validate this on `tutorials/rigidBodyAndBeam`.

## What Was Already Working

Before the latest step, `BlockEigenSolverOF` already solved the rigid-body
translation and rotational correction as extra unknowns. A diagnostic file was
also being written:

```text
tutorials/rigidBodyAndBeam/postProcessing/BlockEigenSolve_rigidBodyMotion/0/BlockEigenSolve_rigidBodyMotion.dat
```

The comparison against:

```text
tutorials/rigidBodyAndBeam/postProcessing/sixDoF_History/0/sixDoFRigidBodyStateFvBeam.dat
```

showed that the BlockEigen rigid-body displacement matched the original sixDoF
solution when the same forcing was used. This gave confidence that the Newmark
time integration terms in `BlockEigenSolverOF` were consistent with the existing
rigid-body solver.

## Latest Implementation

### 1. Extended The Shared Rigid-Body Solution Container

Changed:

```text
src/wireBunchingModels/beamModels/coupledTotalLagNewtonRaphsonBeam/RigidBodyCouplingData.H
```

`RigidBodySolution` now stores more than displacement and rotation correction.
It also stores:

```text
velocity
angularMomentum
acceleration
torque
```

Reason:

The moorFV rigid-body state needs a complete state update, not just position.
If only displacement is copied back, the next time step would still use stale
velocity, acceleration, angular momentum, and torque.

### 2. Stored The Full BlockEigen Rigid-Body Solution

Changed:

```text
src/wireBunchingModels/numerics/BlockEigenSolvers/BlockEigenSolverOF/BlockEigenSolverOF.C
src/wireBunchingModels/numerics/BlockEigenSolvers/BlockEigenSolverOF/BlockEigenSolverOF.H
```

`BlockEigenSolverOF::solve(...)` now copies the computed rigid-body velocity,
angular momentum, acceleration, and torque into `RigidBodySolution`, alongside
the solved displacement and rotation correction.

Reason:

This turns the BlockEigen output into an actual rigid-body update candidate,
rather than just a diagnostic value written to `postProcessing`.

### 3. Exposed The Latest Rigid-Body Solution From beamFoam

Changed:

```text
src/wireBunchingModels/beamModels/beamModel/beamModel.H
src/wireBunchingModels/beamModels/coupledTotalLagNewtonRaphsonBeam/coupledTotalLagNewtonRaphsonBeam.H
src/wireBunchingModels/beamModels/coupledTotalLagNewtonRaphsonBeam/coupledTotalLagNewtonRaphsonBeam.C
src/wireBunchingModels/beamModels/coupledTotalLagNewtonRaphsonBeam/coupledTotalLagNewtonRaphsonBeamEvolve.C
```

The base `beamModel` now has:

```cpp
virtual bool getRigidBodySolution(RigidBodySolution&) const;
```

The `coupledTotalLagNewtonRaphsonBeam` implementation stores the latest
`RigidBodySolution` after the BlockEigen solve and returns it through that
interface.

Reason:

moorFV only sees a generic `beamModel`. It should not need to know the concrete
beam model type just to retrieve the latest rigid-body solution.

### 4. Made The moorFV Restraint Able To Apply The BlockEigen State

Changed in moorFV:

```text
../sDoFRGBFvBeam/sixDoFRigidBodyMotionFvBeam/sixDoFRigidBodyMotionFvBeam.H
../sDoFRGBFvBeam/sixDoFRigidBodyMotionFvBeam/sixDoFRigidBodyMotionFvBeam.C
../sDoFRGBFvBeam/sixDoFRigidBodyMotionFvBeam/restraints/sixDoFRigidBodyMotionFvBeamRestraint/sixDoFRigidBodyMotionFvBeamRestraint.H
../sDoFRGBFvBeam/sixDoFRigidBodyMotionFvBeam/restraints/finiteVolumeBeamRestraint/finiteVolumeBeam.H
../sDoFRGBFvBeam/sixDoFRigidBodyMotionFvBeam/restraints/finiteVolumeBeamRestraint/finiteVolumeBeam.C
```

The restraint API was changed so `restrain(...)` receives a mutable
`sixDoFRigidBodyMotionFvBeam& motion`.

`sixDoFRigidBodyMotionFvBeam` now has:

```cpp
void applyBlockEigenSolution(const RigidBodySolution& rigidBodySolution);
```

This method updates:

```text
centreOfRotation
velocity
acceleration
angularMomentum
torque
orientation
```

The orientation is updated using the existing rigid-body `rotate(...)` logic and
the BlockEigen rotational correction.

Reason:

The rigid-body solver owns the authoritative motion state. Therefore the cleanest
place to apply the BlockEigen state is inside `sixDoFRigidBodyMotionFvBeam`, not
inside beamFoam.

### 5. Added A Switch To Make BlockEigen Authoritative

Changed:

```text
tutorials/rigidBodyAndBeam/constant/dynamicMeshDict
```

Added under the `finiteVolumeBeam` restraint:

```cpp
blockEigenAuthoritativeMotion true;
```

When this switch is true, `finiteVolumeBeam` does the following after
`beam.evolve()`:

1. Retrieves the latest `RigidBodySolution` from the beam model.
2. Calls `motion.applyBlockEigenSolution(...)`.
3. Sets the returned restraint force and moment to zero.

Reason:

The BlockEigen solution already includes the beam reaction contribution. If the
old beam restraint force were also returned to `sixDoFRigidBodyMotionFvBeam`,
the beam force would be applied twice.

### 6. Updated moorFV Include Paths

Changed in moorFV:

```text
../sDoFRGBState/Make/options
```

The `sDoFRGBState` library now includes the beamFoam rigid-body coupling headers.

Reason:

`sixDoFRigidBodyMotionFvBeam.H` now refers to `RigidBodySolution`, which is
declared in beamFoam.

### 7. Case And Diagnostic Scripts

Changed:

```text
tutorials/rigidBodyAndBeam/constant/beam/beamProperties
tutorials/rigidBodyAndBeam/constant/dynamicMeshDict
tutorials/plot_rigidBodyAndBeam_new_vs_old.py
```

The case keeps:

```cpp
blockEigenForceCoupling true;
blockEigenKinematicCoupling false;
```

The comparison script was added to compare the new `rigidBodyAndBeam` case
against `rigidBodyAndBeam_coupledSolver`, including rigid-body displacement and
beam force outputs.

Reason:

This lets the new authoritative BlockEigen path be compared against the previous
coupled-solver algorithm as a sanity check.

## Validation Completed

The following builds passed:

```bash
cd /Volumes/OpenFOAM/colmmcalister-v2312/run/moorFV/src/beamFoam
./Allwmake

cd /Volumes/OpenFOAM/colmmcalister-v2312/run/moorFV
./Allwmake
```

The following case was cleaned and run fresh:

```bash
cd /Volumes/OpenFOAM/colmmcalister-v2312/run/moorFV/src/beamFoam/tutorials/rigidBodyAndBeam
./Allclean
./Allrun
```

The case completed normally to:

```text
End
```

The log confirmed that the new authoritative path was active:

```text
finiteVolumeBeam using authoritative BlockEigen rigid-body solution
```

The regression check also confirmed that the final BlockEigen displacement,
when added to the initial centre of rotation, matches the final sixDoF centre of
rotation within tolerance.

At `t = 0.25`, the final BlockEigen displacement was:

```text
(-0.0001205323727 -5.640956454e-11 -0.2923609842)
```

The final sixDoF centre of rotation was:

```text
(2.6563794678 -5.6383999633e-11 -0.29236098424)
```

This is consistent with the initial centre of rotation:

```text
(2.6565 0 0)
```

because:

```text
2.6565 + (-0.0001205323727) = 2.6563794676273
```

## Current Interpretation

The BlockEigen solution is now authoritative for rigid-body motion in the
`rigidBodyAndBeam` case when `blockEigenAuthoritativeMotion true` is set.

However, this is still not the fully implicit beam-rigid-body coupled matrix.
The current working path uses:

```cpp
blockEigenForceCoupling true;
blockEigenKinematicCoupling false;
```

That means the beam reaction is added into the rigid-body RHS explicitly, and
then the resulting BlockEigen rigid-body solution is written back into the
moorFV rigid-body state.

The important improvement is that the old sixDoF solver is no longer separately
applying the finite-volume beam restraint force when the BlockEigen state is
made authoritative.

## Why This Was Done This Way

The implementation was deliberately switch-gated. This avoids silently changing
the behaviour of other moorFV cases.

The state application was placed inside `sixDoFRigidBodyMotionFvBeam` because
that class owns the rigid-body state and already has the correct orientation
update utilities.

The force returned by `finiteVolumeBeam` is zeroed only when the authoritative
BlockEigen path is active. This preserves the old algorithm when the switch is
disabled.

The beam model exposes the solution through a generic virtual function so moorFV
does not need to cast to `coupledTotalLagNewtonRaphsonBeam`.

## Remaining Risks And Open Questions

1. The current authoritative path is explicit in the beam force contribution.
   It is not yet a fully populated off-diagonal beam-rigid-body block matrix.

2. `blockEigenKinematicCoupling` remains disabled. Earlier tests showed that
   enabling it caused problematic behaviour, likely because the attachment
   boundary motion and matrix kinematic coupling were both trying to control the
   same beam-end degrees of freedom.

3. The log output has been reduced to keep the useful authoritative BlockEigen
   confirmation without the earlier per-iteration vector dumps.

4. The restraint API now takes a mutable `sixDoFRigidBodyMotionFvBeam&`. That is
   needed for this path, but it changes the restraint contract and should be
   reviewed before using this design more broadly.

5. The current test is a single-beam, single-body case. Multi-beam and parallel
   cases still need separate validation.

## Recommended Next Steps

1. Re-run the comparison script.
   Compare `rigidBodyAndBeam` against `rigidBodyAndBeam_coupledSolver` after the
   authoritative switch is enabled.

2. Continue investigating `blockEigenKinematicCoupling`.
   The next physics step is to remove the remaining explicitness by correctly
   populating the off-diagonal matrix terms between the beam-end DOFs and the
   rigid-body DOFs.

3. Decide ownership of attachment kinematics before enabling kinematic coupling.
   The current working state keeps the beam-end boundary prescribed by the
   rigid-body motion. A future matrix-owned mode must avoid prescribing the same
   attachment kinematics twice.

4. Test a longer and more realistic moorFV case.
   After the single-beam tutorial remains stable, move to a multi-beam mooring
   case and verify that the authoritative state update still behaves correctly.
