# `beamFoam` - Geometrically exact Simo–Reissner beam models implemented in OpenFOAM

This repository contains a cell-centred finite volume implementation of 3-D geometrically exact Simo–Reissner beam models in OpenFOAM. 

The default `openfoam-v2306` is tested with OpenFOAM-v2306 but likely works with similar OpenFOAM versions from OpenFOAM.com.


# How to install `beamFoam`

Run `./Allwclean` and `./Allwmake` scripts 


# How to run a test case

To run a test case, execute the `./Allclean` and `./Allrun` scripts located in the case.


# Post-processing using ParaView

The beam's deformed configuration can be observed using the 'Warp by Vector' filter with the `pointW` (point displacement) field.


# Notes on the solver

1. The primary variables are displacement (W) and rotation (Theta) - 6 scalar degrees of freedom per cell.

2. This solver does not use `blockMesh`, but uses its own simple meshing utility, which can be found in `applications/utilities/`.

3. As a starting point, see the test cases in `run/beamFoam/elastic/myCases/AllTestCasesInPaper/`.

Further details can be found in the article [A cell-centered finite volume formulation of geometrically exact Simo–Reissner beams with arbitrary initial curvatures, https://doi.org/10.1002/nme.6994](https://doi.org/10.1002/nme.6994); the test cases from this article are included in the repository.

# Contact

- seevani.bali@ucd.ie
- zeljko.tukovic@fsb.hr
- philip.cardiff@ucd.ie
- amirhossein.taran@ucdconnect.ie
