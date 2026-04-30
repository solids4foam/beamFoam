# Large-Deflection Cantilever Under Self Weight

## Overview
This tutorial reproduces the cantilever beam experiment described by
Belendez, Neipp, and Belendez (2003). The reference problem is a thin,
rectangular steel cantilever that undergoes geometrically nonlinear bending
under its own weight and, optionally, a vertical concentrated load at the free
end.

The default setup in this tutorial is the self-weight bending case. Gravity is
active through `constant/g`, while the free-end displacement boundary condition
in `0/W` applies zero external end force.

## Geometry and Material
- Beam length: `L = 0.4 m`
- Rectangular cross-section width: `b = 0.025 m`
- Rectangular cross-section height: `h = 0.0004 m`
- Discretisation: `80` beam segments
- Young's modulus: `E = 194.3 GPa`
- Shear modulus: `G = 74.73 GPa`
- Density: `rho = 7730 kg/m3`
- Gravity: `(0 0 -9.81) m/s2`

These values are set in `constant/beamProperties` and `constant/g`.

## Boundary Conditions
- Left end: clamped displacement and rotation through `fixedValue` entries in
  `0/W` and `0/Theta`.
- Right end: zero applied moment in `0/Theta`.
- Right end displacement: the default `0/W` setup uses
  `forceBeamDisplacementNR` with `force uniform (0 0 0)`, so the beam bends
  only due to self weight.

To apply a different end force, modify both:
- `constant/timeVsFollowerForce`, which defines the time history and direction
  of the follower force.
- The right-patch boundary condition in `0/W`, enabling or updating the
  `followerForceBeamDisplacementNR` setup so it reads that force series.

If these are not changed, the tutorial should be interpreted as the self-weight
bending case.

## Running the Case
From this directory, run:

```sh
./Allclean
./Allrun
```

`Allrun` creates the beam mesh with `createBeamMesh` and then runs
`beamFoam`.

## Post-Processing
Open the generated case in ParaView and visualise the deformed beam with
`pointW` using `Warp By Vector`. The free-end vertical displacement can be
compared against the reference paper values for the self-weight case or for
cases where a vertical end load is enabled.

## Reference
Belendez, T., Neipp, C., and Belendez, A. (2003). *Numerical and Experimental
Analysis of a Cantilever Beam: a Laboratory Project to Introduce Geometric
Nonlinearity in Mechanics of Materials*. International Journal of Engineering
Education, 19(6), 885-892.
