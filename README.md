# `beamFoam` - Geometrically exact Simo–Reissner beam models implemented in OpenFOAM
This repository contains a cell-centred finite volume implementation of 3-D geometrically exact Simo–Reissner beam models in OpenFOAM.
The default `openfoam-v2306` is tested with OpenFOAM-v2306 but likely works with similar OpenFOAM versions from OpenFOAM.com.


# How to install `beamFoam`
Run `./Allwclean` and `./Allwmake` scripts located in the top-level folder.

# Thirdparty requirements
This repository uses Eigen library third party package. It has been tested with Eigen 3.3.7 release and there is a `Allwmake` script to install Eigen, if not found. the repository likely works with newer releases of Eigen too but it has not been tested.

# How to run a test case
To run a test case, execute the `./Allclean` and `./Allrun` scripts located in the case. Each test case has its own README.md file describing the case.

# Generating plots
Plotting scripts compatible with `gnuplot` utility are located in each test case. Execute the plots using command,
```gnuplot <file-name> ```

# Post-processing using ParaView
The beam's deformed configuration can be observed using the 'Warp by Vector' filter with the `pointW` (point displacement) field.

# Notes on the solver
1. The primary variables are displacement (W) and rotation (Theta) - 6 scalar degrees of freedom per cell.
2. This solver does not use `blockMesh`, but uses its own simple meshing utility `createBeamMesh`, which can be found in `applications/utilities/`.
