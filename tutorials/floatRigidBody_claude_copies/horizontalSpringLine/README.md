# floatRigidBody_claude_copies/horizontalSpringLine

A floating box moored by a single horizontal beamFoam mooring line with a
linear spring at the anchor, in a small wave tank, run in serial with
interFoam and the `sixDoFRigidBodyMotionFvBeam` motion solver.

    ./Allrun      # beams, mesh, setFields, interFoam
    ./Allclean

## Geometry

| Item            | Value                                                      |
|-----------------|------------------------------------------------------------|
| Tank            | 1.0 x 0.05 x 0.3 m, 100 x 5 x 60 cells (10 x 10 x 5 mm)     |
| Water depth     | 0.15 m (`system/setFieldsDict`)                            |
| Box             | 0.1 x 0.05 x 0.05 m, centre (0.5 0.025 0.15), full width   |
| Box draft       | 0.025 m (half submerged)                                   |
| Line anchor     | (0.25 0.025 0.125), linear spring k = 4 N/m along x        |
| Line attachment | box bottom-left corner (0.45 0.025 0.125)                  |

The box is cut out of the blockMesh with `topoSet` + `subsetMesh`; the
exposed faces go into the (initially empty) `floatingObject` wall patch. It
spans the tank width, so the motion is limited to surge, heave and pitch
(plane and axis constraints in `constant/dynamicMeshDict`).

## Mooring line

Set up like the lines in AssessmentOfCoupledFVFramework/caseFiles/baseCase:
a stiff nylon line (R = 1 mm, E = 2.7 GPa, G = 0.97 GPa) whose anchor end
(patch `left`) is a `linearSpringForceBeamDisplacementNR` boundary. The
spring force is -k (W . d) d with d = (1 0 0), so the spring acts only along
the line and the anchor end is free to move sideways. k = 4 N/m is the
Assessment case's 58.7 N/m Froude-scaled by mass (7.12 kg -> 0.1248 kg),
which gives a surge period of roughly 1.2 s.

The line is horizontal, from (0.25 0.025 0.125) to the box's bottom-left
corner (0.45 0.025 0.125), with no pretension (strain 0; with one line
nothing would balance it). It is weightless: as in the Assessment case,
there is no `constant/beamone/g`, so the beam model applies no gravity.
At rest it exerts no force, so the box mass equals its buoyancy.

`beamInit/makeBeams` meshes and places the line (K, E, G are set there). The
beam convergence tolerances `residualTol`, `solutionTol` and `absoluteTol`
are 1e-10, as in the Assessment case: with the default absoluteTol (1e-30)
a line at rest never converges.

## Coupling stability

The box (0.125 kg) is lighter than its heave added mass (about 0.2 kg), so
explicit coupling diverges. The case uses `nOuterCorrectors 3` with
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
