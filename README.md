# beamFoam

beamFoam is a cell-centred finite volume solver for nonlinear geometrically
exact Simo-Reissner beams in OpenFOAM®. It models slender structures undergoing
large three-dimensional deformations and rotations, including tension, bending,
shear and torsion.

beamFoam is a sister repository of
[solids4foam](https://github.com/solids4foam/solids4foam).

## Breaking change: beam boundary condition names

Since [pull request #61][pr61], two beam boundary conditions have new names:

| Old name | New name |
| --- | --- |
| `fixedDisplacement` | `fixedBeamDisplacement` |
| `fixedRotation` | `fixedBeamRotation` |

They were renamed because solids4foam has its own classes with the same names,
and the two clashed when both libraries were loaded ([issue #56][issue56]).
The old names are removed rather than kept as aliases, because an alias would
bring the clash back. Only the names change; the conditions and their entries
are the same.

A case that still uses an old name stops at start-up with an error that names
the patch and the new type. To update a case, run this from its top directory:

```bash
find . -type f \( -name W -o -name Theta \) -exec sed -i \
    -e 's/\bfixedDisplacement\b/fixedBeamDisplacement/g' \
    -e 's/\bfixedRotation\b/fixedBeamRotation/g' {} +
```

This changes only the beam fields `W` and `Theta`, including those in time
directories and `processor*` directories, so restarts keep working. Do not run
a find-and-replace over the whole case: solids4foam's own `fixedDisplacement`
and `fixedRotation`, used on solid fields such as `D`, keep their names.

[pr61]: https://github.com/solids4foam/beamFoam/pull/61
[issue56]: https://github.com/solids4foam/beamFoam/issues/56

## Documentation

The beamFoam user documentation is hosted on the solids4foam website:

- [beamFoam overview](https://solids4foam.github.io/sister-repositories/beamFoam/)
- [installation guide](https://solids4foam.github.io/sister-repositories/beamFoam/installation/)
- [user documentation](https://solids4foam.github.io/sister-repositories/beamFoam/documentation/)
- [tutorial guide](https://solids4foam.github.io/sister-repositories/beamFoam/tutorials/)

## Installation

Source a supported OpenFOAM.com environment, then run the top-level build
scripts:

```bash
./Allwclean
./Allwmake
```

The build script currently recognises OpenFOAM versions `v2106` through `v2506`
at six-month release intervals. OpenFOAM `v2306` is the principal tested
version described in the beamFoam paper.

Eigen is provided under `ThirdParty/` and is built by `Allwmake` when required.

## Running a Tutorial

Each tutorial contains `Allclean` and `Allrun` scripts. For example:

```bash
cd tutorials/cantilever
./Allclean
./Allrun
```

The solver uses the `createBeamMesh` utility rather than `blockMesh`. Cases that
start from a translated or rotated reference configuration also use
`setInitialPositionBeam`.

To visualise a deformed beam in ParaView, apply **Warp By Vector** using the
`pointW` point-displacement field. Many tutorial cases also provide gnuplot
scripts.

## Publication

The formulation, implementation and principal benchmark cases are described in:

> S. Bali, A. Taran, Ž. Tuković, V. Pakrashi and P. Cardiff. beamFoam: A
> Cell-Centred Finite Volume Solver for Nonlinear Geometrically-Exact Beams in
> OpenFOAM. OpenFOAM® Journal, 5, 180-210, 2025.
> [doi:10.51560/ofj.v5.170](https://doi.org/10.51560/ofj.v5.170).

## Support

Please report bugs and request features through the
[beamFoam issue tracker](https://github.com/solids4foam/beamFoam/issues).

## OpenFOAM Trademark

This offering is not approved or endorsed by OpenCFD Limited, producer and
distributor of the OpenFOAM software via the
[OpenFOAM website](https://www.openfoam.com), and owner of the OPENFOAM® and
OpenCFD® trade marks.
