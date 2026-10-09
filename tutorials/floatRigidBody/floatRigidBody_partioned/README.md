# floatRigidBody_partioned

A floating box moored by one beamFoam line in waves, with the box and the
line coupled partitioned (moorFV solver `FvBeamNewmark`: the body and
beamFoam take turns). The line has weight, rests on the tank floor through
seabed contact, and the beam function objects write the line end forces and
displacements. `../floatRigidBody_monolithic` is the same case solved
monolithically; `../floatRigidBody_moorDyn` is the MoorDyn reference.
The same case without seabed contact is
`../trial_cases/floatRigidBody_partioned_gravity_noContact`.

## Gravity on the line

The beam model reads gravity from `constant/<line>/g` and silently uses
g = 0 if that file is missing. `beamInit/makeBeams` did not copy it, so in
the earlier coupled runs the line had weight only during the beamInit
pretension step and was weightless in interFoam. `makeBeams` here also
copies `constant/g` (check: `log.interFoam` contains
"g is set for gravitational body force").

## Shift to the MoorDyn geometry

The case is shifted down 0.15 m to match the MoorDyn case
(foamMooring/floatRigidBody_moorDyn): floor z = -0.15, still water z = 0,
box centre (0.5 0.025 0), attachment (0.45 0.025 -0.025), anchor
(0.25 0.025 -0.15), free-surface probes at z = 0 (blockMeshDict,
topoSetDict, setFieldsDict, controlDict, beamInit/makeBeams,
dynamicMeshDict). The z values in the tables below are those of the
unshifted case. mass (0.124319) balances the line weight at rho = 141.8
(MoorDyn case: 0.124316).

## Running

Run in serial with interFoam and the `sixDoFRigidBodyMotionFvBeam` motion
solver (about 65 min on one core):

    ./Allrun      # beams, mesh, setFields, interFoam
    ./Allclean

## Seabed contact

Seabed contact is switched on only by `constant/beamMomentumContributionProperties`
(type `groundContact`, read by every beam in the case). The
`groundContactActive`, `groundZ`, `gStiffness` and `gDamping` keys in
`beamProperties` are not read by the code.

| Setting      | Value                                   |
|--------------|-----------------------------------------|
| groundZ      | -0.15 m (tank floor, MoorDyn `WtrDpth`) |
| kNormal      | 1e4 Pa/m                                |
| cNormal      | 1 Pa s/m                                |
| Friction     | none (as in MoorDyn)                    |

Contact acts on the cell centres below `groundZ`, with a force per unit
length max(2R (kNormal penetration - cNormal Uz), 0). The normal stiffness
is in the Jacobian; damping is explicit. kNormal = 1e4 was chosen before the
contact force was fixed (8 Oct 2026: wrong sign, no cell length, no
Jacobian, so the first floor contact diverged); MoorDyn's kBot = 3e6 has not
been tried since.

Check in `log.interFoam`: "Found beamMomentumContribution type:
groundContact" once, and "Number of cells in contact : N" in every beam
iteration. The line first touches the floor at t = 2.60 s, after the snap
load, and up to 15 cells are then in contact. The floor carries the slack
line: the mean axial force from t = 6 s is -0.54 mN at the anchor and
+0.004 mN at the box (without contact: -0.22 mN and +0.33 mN). Peak tension
at the box: 0.1006 N at t = 2.316 s.

## Output

As well as moorFV's restraint output in `postProcessing/0`
(`axialForcebeamone.dat`: time, anchor axial, anchor shear, attachment axial,
attachment shear; `anchorForcebeamone.dat`, `attachmentForcebeamone.dat`),
four beam function objects run from the main `system/controlDict` with
`region beamone;` (in a coupled run `system/beamone/controlDict` is not
read):

| Function object             | File in postProcessing/0                  |
|-----------------------------|--------------------------------------------|
| beamForcesMomentsAnchor     | beamForcesMoments_left.dat (anchor)        |
| beamForcesMomentsAttachment | beamForcesMoments_right.dat (box)          |
| beamDisplacementsAttachment | beamDisplacements_right.dat (box)          |
| beamConvergenceData         | beamConvergenceData.dat                    |

The force Q in beamForcesMoments and in anchor/attachmentForce is the end
force vector (axial plus shear), not the tension: once the line is slack
and bends, most of it is shear. The tension is the axial force in
`axialForcebeamone.dat`. With outer correctors these files have several
rows per time step; keep the last one.

Plot with `../../plot_monolithic_vs_partitioned_floatRigidBody.py`,
`../../plot_floatRigidBody_moorDyn_comparison.py`,
`../../plot_monolithic_vs_partitioned_motion.py` and
`../../plot_monolithic_vs_partitioned_cost.py`.

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

## Waves

StokesII, H = 0.02 m, T = 1.0 s, ramped over 2 s (`constant/waveProperties`),
shallow-water absorption at the outlet. The flow is laminar. Free-surface
probes at x = 0.25 and 0.75 m (`postProcessing/interfaceHeight1`).
