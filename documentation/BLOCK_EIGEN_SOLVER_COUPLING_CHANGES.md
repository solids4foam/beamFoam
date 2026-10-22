# BlockEigenSolverOF Coupling Changes

This note explains the recent changes made to
`src/wireBunchingModels/numerics/BlockEigenSolvers/BlockEigenSolverOF/BlockEigenSolverOF.C`
and the related wiring in `coupledTotalLagNewtonRaphsonBeamEvolve.C`.

## Background

`BlockEigenSolverOF` converts the beamFoam 6-DOF block matrix into an Eigen sparse
matrix. The matrix already contained the beam unknowns:

- `W_x, W_y, W_z`
- `Theta_x, Theta_y, Theta_z`

The coupled rigid-body work appends six additional unknowns:

- rigid-body displacement correction, 3 components
- rigid-body rotation correction, 3 components

Before these changes, the appended rigid-body block was mostly an identity block.
The beam rows could depend on rigid-body unknowns through kinematic coupling, but
the rigid-body rows did not contain reciprocal beam-to-rigid-body terms.

In matrix form, the old structure was approximately:

```text
[ K_beam      K_beam-rb ] [ beam unknowns ] = [ beam RHS ]
[ 0           I         ] [ rb unknowns   ] = [ rb RHS   ]
```

This meant the rigid body could receive beam forces only through explicit RHS
updates. Those beam forces were not represented as implicit matrix coupling.

## What Changed

The appended rigid-body block remains an identity block. This is intentional:
the appended rigid-body unknowns are displacement and rotation updates, not
accelerations. The Newmark-beta displacement/rotation updates are therefore
placed directly in the rigid-body RHS.

The translational RHS uses:

```text
x_(n+1) = x_n + v_n*deltaT
        + deltaT^2*(beta*a_(n+1) + (0.5 - beta)*a_n)
```

The rotational RHS follows the same current code structure, using the previous
angular momentum and current/previous torque terms:

```text
thetaCorrection = deltaT*pi_n
                + deltaT^2*(beta*tau_(n+1) + (0.5 - beta)*tau_n)
```

The matrix diagonal for these six rigid-body rows is therefore:

```text
I_6x6
```

not `mass/(beta*deltaT^2)` or `I/(beta*deltaT^2)`.

The matrix now also receives reciprocal beam-to-rigid-body coupling terms. The
structure is closer to:

```text
[ K_beam      K_beam-rb ] [ beam unknowns ] = [ beam RHS ]
[ K_rb-beam   I         ] [ rb unknowns   ] = [ rb RHS   ]
```

The new rigid-body force-row coefficients use the same local beam boundary
Jacobian already used by the beam attachment coupling:

```cpp
rigidBodyBeamForceWCoeff = Cw/pDelta;
rigidBodyBeamForceThetaCoeff = CQTheta_.boundaryField()[patchI][faceI];
```

where:

- `Cw` is the boundary `CQW_` tensor for force sensitivity to displacement.
- `pDelta` is the boundary-cell-to-face distance, computed from
  `1.0/mesh().deltaCoeffs()`.
- `CQTheta_` is the force sensitivity to beam rotation.

The moment arm is now:

```cpp
momentArm = localRigidBodyData.attachmentPoint
          - localRigidBodyData.centreOfRotation;
```

This is the vector from the rigid-body centre of rotation to the beam attachment
point. It is used to add the rotational rigid-body row coupling from the force
Jacobian through:

```text
torque contribution = r x force contribution
```

In code, this is assembled with a local spin tensor helper:

```cpp
localSpinTensor(momentArm) & forceJacobian
```

The rigid-body rows therefore now receive matrix terms from:

- beam displacement unknowns to rigid-body force rows
- beam rotation unknowns to rigid-body force rows
- beam displacement unknowns to rigid-body torque rows through `r x F`
- beam rotation unknowns to rigid-body torque rows through `r x F`

## Why These Changes Were Made

The previous implementation was only partly coupled in the matrix. It allowed
the beam equation to depend on rigid-body motion, but the rigid-body equation did
not implicitly depend on the beam solution.

That is not a fully coupled beam/rigid-body linear system. It is closer to a
one-way matrix coupling plus explicit force feedback.

The goal of these changes was to move the beam-to-rigid-body force feedback into
the same linear system. This should make the Eigen solve represent the coupled
beam and rigid-body response more directly, instead of relying only on explicit
RHS force updates.

## Files Touched

- `BlockEigenSolverOF.H`
  Added storage and constructor arguments for the reciprocal beam force
  Jacobians and the rigid-body moment arm.

- `BlockEigenSolverOF.C`
  Keeps the rigid-body block as identity, puts the Newmark-beta rigid-body
  displacement/rotation updates in the RHS, and adds reciprocal
  beam-to-rigid-body matrix entries plus a local spin tensor helper.

- `coupledTotalLagNewtonRaphsonBeamEvolve.C`
  Computes and passes `Cw/pDelta`, `CQTheta_`, and the correct current
  `momentArm` into `BlockEigenSolverOF`.

## Verification

The code was rebuilt successfully with:

```bash
./Allwmake
```

from `src/beamFoam`, and then the top-level moorFV libraries were rebuilt with:

```bash
BEAMFOAM_DIR=/Users/colmmcalister/OpenFOAM/colmmcalister-v2312/run/moorFV/src/beamFoam ./Allwmake
```

Both builds completed successfully. The build still emits existing Eigen/header
warnings, but no compile errors.

## Remaining Caveat

This has been compile-verified only. A coupled tutorial or benchmark case should
be run next to check the numerical behavior, especially the signs of the new
reciprocal force and torque Jacobian terms.
