# floatRigidBody_partioned

A floating box moored by a single beamFoam mooring line in a small wave tank,
run in serial with interFoam and the `sixDoFRigidBodyMotionFvBeam` motion
solver.

    ./Allrun      # beams, mesh, setFields, interFoam
    ./Allclean

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

The box (0.09 kg) is lighter than its heave added mass (about 0.2 kg), so
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

## Waves

StokesII, H = 0.02 m, T = 1.0 s, ramped over 2 s (`constant/waveProperties`),
shallow-water absorption at the outlet. The flow is laminar. Free-surface
probes at x = 0.25 and 0.75 m (`postProcessing/interfaceHeight1`).
