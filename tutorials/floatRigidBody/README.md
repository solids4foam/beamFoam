# floatRigidBody

A floating box moored by two beamFoam mooring lines in a small wave tank,
run in serial with interFoam and the `sixDoFRigidBodyMotionFvBeam` motion
solver. `../floatRigidBody_claude` is the same case with a single,
unpretensioned line.

    ./Allrun      # beams, mesh, setFields, interFoam
    ./Allclean

## Geometry

| Item            | Value                                                      |
|-----------------|------------------------------------------------------------|
| Tank            | 1.0 x 0.05 x 0.3 m, 100 x 5 x 60 cells (10 x 10 x 5 mm)     |
| Water depth     | 0.15 m (`system/setFieldsDict`)                            |
| Box             | 0.1 x 0.05 x 0.05 m, centre (0.5 0.025 0.15), full width   |
| Box draft       | 0.025 m (half submerged)                                   |
| Line anchors    | (0.25 0.025 0) and (0.75 0.025 0)                          |
| Line attachment | box bottom corners (0.45 0.025 0.125), (0.55 0.025 0.125)  |

The box is cut out of the blockMesh with `topoSet` + `subsetMesh`; the
exposed faces go into the (initially empty) `floatingObject` wall patch. It
spans the tank width, so the motion is limited to surge, heave and pitch
(plane and axis constraints in `constant/dynamicMeshDict`).

## Mooring lines

`beamInit/makeBeams` builds each line from `beamInit/template`:

1. `createCircularBeamMesh` at the unstretched length L0 = L/(1 + strain).
2. `transformPoints` + `setInitialPositionBeam` to place it anchor -> box.
3. `beamFoam` (steady) pulls the attachment end out to the box corner.
4. The final state is copied to `0.orig/<beam>`, with `constant/<beam>` and
   `system/<beam>`; W on the attachment patch becomes fixedValue (driven by
   the rigid body in the coupled run) and steadyState is switched off.

Line properties: R = 1 mm, E = 5 MPa, rho = 1140 kg/m3, strain 0.02, giving
about 0.315 N tension per line (0.1707 N vertical at the box).

The box mass is buoyancy minus the vertical pull of the lines, so the box
starts in equilibrium. **If you change the line geometry, E, R or strain,
recompute `mass` in `constant/dynamicMeshDict`** from the new vertical
attachment force (z-component of `Q` on the `right` patch in
`0.orig/<beam>/Q`).

Note: the beam model currently applies the line's dry weight (the buoyancy
term in coupledTotalLagNewtonRaphsonBeamEvolve.C is commented out).

## Coupling stability

The box (0.09 kg) is lighter than its heave added mass (about 0.2 kg), so
explicit coupling diverges. The case uses `nOuterCorrectors 3` with
`moveMeshOuterCorrectors yes` and `accelerationRelaxation 0.4`. With these
settings the case ran stably to 6.15 s in waves (with kEpsilon, before the
switch to laminar; the laminar version has not been run to completion).

## Restarting

The beam regions do not write their reference geometry (`refW`, `refWf`,
`refLambda`, `refLambdaf`, `refTangent`) at output times, and without it the
lines diverge on restart. Copy them in from 0 before restarting:

    for b in beamone beamtwo; do cp 0/$b/ref* <time>/$b/; done

Line velocity and acceleration are not written either, but the effect on a
restart was below 0.2 % of the line force.

## Waves

StokesII, H = 0.02 m, T = 1.0 s, ramped over 2 s (`constant/waveProperties`),
shallow-water absorption at the outlet. The flow is laminar. Free-surface
probes at x = 0.25 and 0.75 m (`postProcessing/interfaceHeight1`).
