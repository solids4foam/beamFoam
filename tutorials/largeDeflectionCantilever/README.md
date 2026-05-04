# Large-Deflection of Cantilever (rectangular cross-section) Under Self Weight

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
- Weight of the beam is `0.3032 N`
- Discretisation: `80 cells`
These values are set in `constant/beamProperties` and `constant/g`.

Two options to add self weight of beam:
1. Use q field  below - a uniform distributed load representing self-weight
or 2. Use density and gravity - can be set in constant/

## For a  static load case:
- Set the self-weight of the beam as uniform distributed load (UDL) of `0.758 N/m`.
  This UDL can be defined by `0/q` volVectorField

## For running a dynamic case, set:
- Density: `rho = 7730 kg/m^3`
- Gravity: `(0 0 -9.81) m/s^2` in `constant/g`
- Set `ddtScheme` and `d2dt2Scheme` in `system/fvSchemes` to Euler or Newmark.
- NOTE: Euler introduces numerical damping unless time-step size is small (e.g.`0.0001`).
  Check the beam energy via the function-objects (see below).
- IMP NOTE: Remove the `q` field else it will take the self-weight into account twice.

## Boundary Conditions
- Left end: clamped displacement and rotation through `fixedValue` entries in
  `0/W` and `0/Theta`.
- Right end: zero applied moment in `0/Theta`.
- Right end displacement: the default `0/W` setup uses
  `forceBeamDisplacementNR` with `force uniform (0 0 0)`, so the beam bends
  only due to self weight.

To apply a different end force, modify both:
- `constant/timeVsForce`, which defines the time history and direction
  of the force.
- The right-patch boundary condition in `0/W`, enabling or updating the
  `forceBeamDisplacementNR` setup so it reads that force series.

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

## Function Objects
This case enables two function objects in `system/controlDict`:

- `beamDisplacements1`
  Writes beam displacement history for patch `right` to
  `postProcessing/0/beamDisplacements_right.dat`.
- `beamEnergyData`
  Writes energy history to `postProcessing/0/beamEnergyData.dat`.

These outputs are useful for dynamic runs:

- `beamDisplacements_right.dat` can be used to track free-end deflection with
  time.
- `beamEnergyData.dat` can be used to monitor internal, kinetic, and total
  energy, which is useful when checking time-integration behaviour. This
  energy output is only valid for dynamic cases.

The provided `allPlots.gnuplot` script reads these two files directly and
generates `energyPlot.pdf` and `displacementPlot.pdf`.

## Post-Processing
Open the generated case in ParaView and visualise the deformed beam with
`pointW` using `Warp By Vector`. The free-end vertical displacement can be
compared against the reference paper values for the self-weight case or for
cases where a vertical end load is enabled.

For quick history plots from the function-object output, run:

```sh
gnuplot allPlots.gnuplot
```

## Reference paper
Belendez, T., Neipp, C., and Belendez, A. (2003). *Numerical and Experimental
Analysis of a Cantilever Beam: a Laboratory Project to Introduce Geometric
Nonlinearity in Mechanics of Materials*. International Journal of Engineering
Education, 19(6), 885-892.

## Validation Against Reference Paper

The numerical results can be directly compared with experimental data from  
Belendez et al. (2003). The beam is discretised with 80 beam cells.

### Table 1: Free-end vertical displacement \( \delta_y \) vs applied load \( F \)

| F (N) | δy (m) Experimental | δy (m) Numerical (E = 194.3 GPa) | Relative Error |
|------:|--------------------:|----------------------------------:|---------------:|
| 0.000 | 0.089  | 0.0898 | 0.89% |
| 0.098 | 0.149  | 0.1516 | 1.71% |
| 0.196 | 0.195  | 0.1960 | 0.25% |
| 0.294 | 0.227  | 0.2270 | 0.51% |
| 0.392 | 0.251  | 0.2495 | 0.60% |
| 0.490 | 0.268  | 0.2659 | 0.78% |
| 0.588 | 0.281  | 0.2784 | 0.90% |

