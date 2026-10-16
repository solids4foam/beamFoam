# floatRigidBody_twoLine/floatRigidBody_partitioned

The floating box of `../../floatRigidBody/floatRigidBody_partioned` with a
second, mirror-image mooring line on the other side of the box. It is the
Stage 0 reference case of `documentation/phase5MultiLinePlan.md`: the
partitioned (moorFV solver `FvBeamNewmark`) two-line solution that the
multi-line monolithic solver will be compared against. `beamFoamCoupled`
cannot run it yet; it accepts only one `finiteVolumeBeam` restraint.

Built 13 October 2026; not yet run.

## Changes from the one-line case

| File                        | Change                                                     |
|-----------------------------|------------------------------------------------------------|
| `beamInit/makeBeams`        | Second row in `lines`: `beamtwo`                           |
| `constant/dynamicMeshDict`  | Restraint `beamtwo`; `mass` 0.124319 -> 0.123863           |
| `system/controlDict`        | Comment only: beam function objects stay on beamone        |
| `Allclean`                  | Removes the generated beamtwo folders too                  |

Everything else (tank, waves, box, line properties, seabed contact, PIMPLE
settings) is unchanged.

## Lines

| Line    | Anchor               | Attachment (box edge)        |
|---------|----------------------|------------------------------|
| beamone | (0.25 0.025 -0.15)   | (0.45 0.025 -0.025), lower left  |
| beamtwo | (0.75 0.025 -0.15)   | (0.55 0.025 -0.025), lower right |

The lines are mirror images about the box centre x = 0.5. Both are built the
same way by `beamInit/makeBeams` (R = 1 mm, E = 5 MPa, rho = 141.8 kg/m3
submerged, 40 segments, strain 0). In each line patch `left` is the anchor
and `right` the attachment to the box, so both restraints use
`attachmentPatch right; anchorPatch left;`.

The strain is kept at 0 so each line matches the one-line case. With two
opposed lines a pretension could now be balanced, but that would change the
line stiffness and the comparison.

## Box mass

Each line's sagging weight pulls down on the box with 0.00447 N
(0.000456 kg * g, from the one-line case), so

    mass = buoyancy - 2 * 0.000456 = 0.124775 - 0.000912 = 0.123863 kg

After `beamInit/makeBeams`, check this against the z-component of `Q` on the
`right` patch of `0.orig/beamone/Q` and `0.orig/beamtwo/Q`, and correct
`mass` in `constant/dynamicMeshDict` if they differ.

The horizontal pulls of the two lines (about 0.0063 N each) are opposed and
cancel, so the slow initial drift of the one-line case should not appear.
The moment of inertia is unchanged (the box is the same).

## Output

The restraint writes one set of files per line in `postProcessing/0`:
`axialForcebeamone.dat` and `axialForcebeamtwo.dat` (time, anchor axial,
anchor shear, attachment axial, attachment shear), and likewise
`anchorForce*`, `attachmentForce*`, `force*` and `displacement*`.

The beam function objects (`beamForcesMoments`, `beamDisplacements`,
`beamConvergenceData`) run only for beamone. Their file names carry the patch
name but not the region, so a beamtwo copy would overwrite beamone's files.
Giving them per-region file names is a small code change in
`src/wireBunchingModels/functionObjects`, needed before Stage 4/5 of the plan.

## Running

    ./Allrun      # beams, mesh, setFields, interFoam (serial)
    ./Allclean

Expect somewhat more than the one-line case's ~65 min on one core, because
there are two line solves per outer corrector. The plan's reference runs use
`nOuterCorrectors` 3 (as set in `system/fvSolution`) and 8.

Check in `log.interFoam` that each line reads gravity ("g is set for
gravitational body force", once per line), that `groundContact` is found, and
that both `beamone` and `beamtwo` restraints are constructed.

## Restarting

As in the one-line case, copy the reference fields into the restart time for
both lines:

    cp 0/beamone/ref* <time>/beamone/
    cp 0/beamtwo/ref* <time>/beamtwo/

See `../../floatRigidBody/floatRigidBody_partioned/README.md` for the
seabed contact, gravity, coupling-stability and wave details shared with
this case.
