# floatRigidBody_moorDyn

The `floatRigidBody_base` case (floating box, one mooring line, small wave
tank) with the beamFoam line replaced by a MoorDyn v2 line, using the
standard `sixDoFRigidBodyMotion` solver and the foamMooring `moorDynR2`
restraint (`libsixDoFMooring`).

    ./Allrun      # mesh, setFields, interFoam
    ./Allclean

## Differences from floatRigidBody_base

| Item | base (beamFoam) | this case (MoorDyn v2) |
|---|---|---|
| Motion solver | `sixDoFRigidBodyMotionFvBeam` | `sixDoFRigidBodyMotion` motion with modified morphing: `sixDoFRigidBodyMotionModifiedMorphing` (+ `libs (sixDoFMooring)`) |
| Mooring | `finiteVolumeBeam` restraint, `beamInit/makeBeams` | `moorDynR2`, POINT coupling, `Mooring/lines.txt` |
| Vertical datum | floor z = 0, still water z = 0.15 | floor z = -0.15, still water z = 0 |
| Box centre | (0.5 0.025 0.15) | (0.5 0.025 0) |
| Fairlead | (0.45 0.025 0.125) | (0.45 0.025 -0.025) |
| Anchor | (0.25 0.025 0) | (0.25 0.025 -0.15) |
| State function object | `sixDoFRigidBodyStateFvBeam` | `sixDoFRigidBodyState` |
| Mesh morphing | moorFV modified morphing (`xDistance`, `yDistance`) | the same, via `sixDoFRigidBodyMotionModifiedMorphing` (moorFV `src/sixDoFModifiedMorphing`) |

Everything is shifted down 0.15 m because MoorDyn puts the still-water level
at z = 0 and the seabed at z = -WtrDpth, and foamMooring passes OpenFOAM
coordinates to MoorDyn unchanged. The wave boundary conditions measure
height from the bottom of the inlet patch, so they are unaffected.

## Mooring line (Mooring/lines.txt)

Same nylon line as the beam case: R = 1 mm, E = 5 MPa, rho = 1140 kg/m3,
unstretched length 0.23585 m (straight anchor-to-fairlead distance, no
pretension), 40 segments.

- Mass/m = rho*pi*R^2 = 3.5814e-3 kg/m. MoorDyn applies buoyancy itself, so
  the real density is used (the beam case used the submerged density).
- EA = E*pi*R^2 = 15.708 N, EI = E*pi*R^4/4 = 3.927e-6 N m2.
- Internal damping 80 % of critical (BA = -0.8), Cd = 1.2, Ca = 1.0,
  CdAx = 0.05.
- dtM = 2e-5 s (RK2); seabed kBot = 3e6 Pa/m, cBot = 3e5 Pa s/m.

MoorDyn finds the initial line shape itself (dynamic relaxation). The
resulting force on the box at t = 0 is (-0.00636 0 -0.00451) N, matching the
beam line (-0.0063 0 -0.00447 N). **If you change the line, recompute `mass`
in `constant/dynamicMeshDict`** from the z-component of the
`t = 0  attachPt[0]  force` line in `log.interFoam`.

The damping, drag and added-mass coefficients are not taken from the beam
case (it has no equivalents); adjust them if you are comparing the two.

## Output

- `Mooring/lines.out`: fairlead and anchor tension (FairTen1, AnchTen1)
- `Mooring/lines_Line1.out`: node positions and segment tensions
- `Mooring/VTK/`: line geometry for ParaView. The `mdv2_NNNN.vtk` files carry
  no time, so `./makeMooringVTKSeries.py` (also run by `Allrun`, and safe to
  run while the case is running) writes `mdv2.vtk.series` with each file's
  time. Open `case.foam` and `Mooring/VTK/mdv2.vtk.series` together. The line
  lies at y = 0.025, inside the tank, so make the case surface transparent
  or slice it to see the line; a Tube filter makes it easier to see.
- `postProcessing/sixDoF_History`, `postProcessing/interfaceHeight1`

## Restarting

The restraint writes `Mooring/restartFile_<time>` at each write time. To
restart, add `restartFile "Mooring/restartFile_<time>";` to the `moorDynR2`
restraint in `constant/dynamicMeshDict`; otherwise MoorDyn restarts from its
initial shape.

## Waves and coupling

Unchanged from the base case: StokesII, H = 0.02 m, T = 1.0 s, ramped over 2 s;
`nOuterCorrectors 3`, `moveMeshOuterCorrectors yes`,
`accelerationRelaxation 0.4`.

## Mesh morphing

The standard sixDoFRigidBodyMotion morphing blends all of the motion over
innerDistance/outerDistance (0.01/0.1 m). With that, the box's drift toward
the anchor (about 5.6 cm by t = 3.75 s) was squeezed into the 9 cm band to
its left: cells there fell to about 1 % of their volume (non-orthogonality
85 deg, aspect ratio 138), the Courant number jumped to 92 and the run
crashed at t = 3.75 s (the floating-point exception surfaced in MoorDyn's
output writer, but the flow had already blown up).

The case now uses `sixDoFRigidBodyMotionModifiedMorphing` from
`moorFV/src/sixDoFModifiedMorphing` (build with `wmake libso` there): the
standard solver with moorFV's modified morphing, as in the beamFoam cases.
The surge is applied rigidly over the box's x range and blended to zero over
`xDistance` (0.3 m) through the whole depth. Its `motionScale` and
`xmotionScale` fields are identical to those of the beamFoam
floatRigidBody_partioned case.
