# floatRigidBody_monolithic_gravity

floatRigidBody_monolithic with gravity on the mooring line in the coupled run.

The beam model reads gravity from `constant/<line>/g` and silently uses
g = 0 if that file is missing. `beamInit/makeBeams` did not copy it, so in
the earlier coupled runs the line had weight only during the beamInit
pretension step and was weightless in interFoam. `makeBeams` here also
copies `constant/g` (check: `log.interFoam` contains
"g is set for gravitational body force").

The case is also shifted down 0.15 m to match the MoorDyn case
(foamMooring/floatRigidBody_moorDyn): floor z = -0.15, still water z = 0,
box centre (0.5 0.025 0), attachment (0.45 0.025 -0.025), anchor
(0.25 0.025 -0.15), free-surface probes at z = 0 (blockMeshDict,
topoSetDict, setFieldsDict, controlDict, beamInit/makeBeams,
dynamicMeshDict). The z values in the tables below are those of the
unshifted case. mass (0.124319) balances the line weight at rho = 141.8
(MoorDyn case: 0.124316).

floatRigidBody_partioned with the box and its mooring line solved monolithically:
moorFV solver `beamFoamCoupled` instead of the partitioned `FvBeamNewmark`.
Everything else (mesh, line, waves, `nOuterCorrectors 3`) is the same, so the
two cases can be compared with
`../../plot_monolithic_vs_partitioned_floatRigidBody.py` (motion, line forces,
free surface) and `../../plot_monolithic_vs_partitioned_cost.py`
(cost).

    ./Allrun      # beams, mesh, setFields, interFoam
    ./Allclean

Differences from floatRigidBody_partioned (`constant/dynamicMeshDict`):

- `solver { type beamFoamCoupled; }`. The plane/axis constraints are the
  same; beamFoamCoupled applies them inside the monolithic solve (below).
- `accelerationRelaxation` is not used: `beamFoamCoupled` accepts beamFoam's
  body state as it is.

Here, "monolithic" describes the box--mooring-line subsystem. The fluid is
still coupled to that subsystem through the three PIMPLE outer correctors.

## Completed-run comparison

The current `floatRigidBody_monolithic` and `floatRigidBody_partioned` results
both finish stably at 10 s with no beam Newton failures. Their authored
numerical inputs are identical apart from the rigid-body/beam solver and the
use of acceleration relaxation by the partitioned solver.

| Quantity at 10 s | Monolithic | Partitioned | Partitioned - monolithic |
|------------------|------------|-------------|----------------------------|
| Surge            | -73.303 mm | -86.344 mm  | -13.041 mm                 |
| Heave            | +5.632 mm  | +6.352 mm   | +0.720 mm                  |
| Pitch            | +2.908 deg | +2.793 deg  | -0.115 deg                 |
| Attachment tension | 0.7550 mN | 0.7508 mN | -0.0042 mN                 |

Over the full histories, the largest differences are 14.09 mm in surge,
0.868 mm in heave and 0.489 deg in pitch. Sway, roll and yaw remain exactly
zero in both cases. The maximum free-surface differences at the two probes
are only 0.317 mm and 0.130 mm, so the main difference is accumulated
horizontal drift rather than the wave field.

The peak attachment tension is 0.0921 N at 1.230 s for the monolithic case
and 0.0964 N at 1.206 s for the partitioned case. The peak magnitudes differ
by 4.6%, with a small phase shift during the transient.

| Cost metric | Monolithic | Partitioned |
|-------------|------------|-------------|
| Time steps | 5001 | 5012 |
| Mean summed beam Newton iterations per time step | 12.758 | 11.827 |
| Execution time | 3214.52 s | 3189.63 s |
| Wall time | 3234 s | 3191 s |

The monolithic case is 0.78% slower by OpenFOAM execution time and 1.35%
slower by wall time in these runs. At three outer correctors the schemes are
therefore similar in cost but are not coupling-converged to the same surge
trajectory. With 8 outer correctors (the `_nOuter8` cases) the final surge
difference falls to 1.3 mm, and the 3 corrector monolithic run is much closer
to the converged drift than the 3 corrector partitioned run. The two methods
nevertheless give completely different mooring-line shapes. See
`../compare_partioned_monolithic_claude.txt` for the full comparison.

## Constraints in the monolithic solve

The constraints are needed, not only convenient: the half-submerged
0.05 x 0.05 m cross-section has a negative roll metacentric height
(GM = KB + BM - KG = 0.0125 + 0.0083 - 0.025 = -0.0042 m), so the box is
hydrostatically unstable in roll. Without constraints, roll grew from about
1e-7 rad at about 15 /s, took sway with it, and the run diverged at
t = 1.07 s. Pitch is stable (GM = +0.021 m).

beamFoamCoupled passes the motion's projections onto the free directions to
beamFoam every step: `tConstraints()` (global axes) and
`Q & rConstraints() & Q.T` (rConstraints is held in body axes). In the
rigidBodyEnd rows and columns of the BlockEigen system
(`constrainRigidBodyEndCoupling`), with P the projection and C = I - P:

- body columns are multiplied by P on the right, body rows by P on the left;
- C, scaled like the block, is added to the body diagonal blocks;
- the body sources (minus the residuals) are projected, P b.

So the constrained increments solve C dx = 0, and the free ones keep their
Newton equations. The constrained parts of the fluid force and moment drop
out of the residual, as they do in FvBeamNewmark. rigidBodyEnd also projects
its state and the partitioned Newmark accelerations, for beamFoam's own
partitioned coupling. The Jacobian check (`jacobianCheck`) still compares
the unprojected blocks.

This was the first monolithic case whose line is placed by
`setInitialPositionBeam`, i.e. with a non-zero reference displacement `refWf`.
rigidBodyEnd used to take the attachment point as `Cf + W` without `refWf`,
which put the arm 0.3 m off and gave spurious roll and yaw; it now uses
`Cf + refWf + W` (coupledTotalLagNewtonRaphsonBeamRigidBodyEnd.C).

## Geometry

| Item            | Value                                                      |
|-----------------|------------------------------------------------------------|
| Tank            | 1.0 x 0.05 x 0.3 m, 100 x 5 x 60 cells (10 x 10 x 5 mm)     |
| Water depth     | 0.15 m (`system/setFieldsDict`)                            |
| Box             | 0.1 x 0.05 x 0.05 m, centre (0.5 0.025 0.15), full width   |
| Box draft       | 0.025 m (half submerged)                                   |
| Line anchor     | (0.25 0.025 0)                                             |
| Line attachment | box bottom-left corner (0.45 0.025 0.125)                  |

The box is cut out of the blockMesh with `topoSet` + `subsetMesh`; the
exposed faces go into the (initially empty) `floatingObject` wall patch. It
spans the tank width, so the motion is limited to surge, heave and pitch
(plane and axis constraints in `constant/dynamicMeshDict`).

## Mooring line

`beamInit/makeBeams` builds the line from `beamInit/template` (add rows to
`lines` in that script for more lines):

1. `createCircularBeamMesh` at the unstretched length L0 = L/(1 + strain).
2. `transformPoints` + `setInitialPositionBeam` to place it anchor -> box.
3. `beamFoam` (steady) pulls the attachment end out to the box corner and
   lets the line settle under its weight.
4. The final state is copied to `0.orig/beamone`, with `constant/beamone`
   and `system/beamone`; W on the attachment patch becomes fixedValue
   (driven by the rigid body in the coupled run) and steadyState is
   switched off.

Line properties: R = 1 mm, E = 5 MPa, strain 0 (no pretension). With one
line nothing would balance a pretension, so the line starts slack-straight
and only takes load when the box moves away from the anchor. It sags about
3 mm under its own weight, giving a tension of about 7.7 mN.

The beam model does not apply buoyancy to the line (the term in
coupledTotalLagNewtonRaphsonBeamEvolve.C is commented out), so `rho` is set
to the submerged density of nylon, 1140 - 998.2 = 141.8 kg/m3.

The box mass is buoyancy minus the vertical pull of the line (0.00447 N),
so the box starts in vertical equilibrium. The line's horizontal pull of
0.0063 N toward the anchor is not balanced, so expect a slow initial drift.
**If you change the line geometry, E, R, rho or strain, recompute `mass` in
`constant/dynamicMeshDict`** from the z-component of `Q` on the `right`
patch in `0.orig/beamone/Q`.

## Coupling stability

The box (0.124319 kg) is lighter than its heave added mass (about 0.2 kg), so
explicit coupling diverges (it was lighter still with the earlier
two-line pretension). The case uses `nOuterCorrectors 3` with
`moveMeshOuterCorrectors yes` and `accelerationRelaxation 0.4`.

## Restarting

The beam regions do not write their reference geometry (`refW`, `refWf`,
`refLambda`, `refLambdaf`, `refTangent`) at output times, and without it the
lines diverge on restart. Copy them in from 0 before restarting:

    cp 0/beamone/ref* <time>/beamone/

Line velocity and acceleration are not written either, but the effect on a
restart was below 0.2 % of the line force.

This restart procedure applies to the partitioned solver. `beamFoamCoupled`
currently requires the initial body state and rejects a restart from a
displaced centre of mass.

## Waves

StokesII, H = 0.02 m, T = 1.0 s, ramped over 2 s (`constant/waveProperties`),
shallow-water absorption at the outlet. The flow is laminar. Free-surface
probes at x = 0.25 and 0.75 m (`postProcessing/interfaceHeight1`).
