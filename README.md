# bali-beamSolver - Seevani's version of the beamsolver code

FVM (OpenFOAM) (cell-centered) implementation of 3-D Simo Reissner Beam Solver

Currently works only with foam-extend-4.1

This branch has many experimental test cases which might not run.

Contact: seevani.bali@ucdconnect.ie if some test case is not working.

# First Compile:

Run ./Allwclean and ./Allwmake scripts 

# Run a test case:

As of now, you can run all test cases under run/beamFoam/elastic/ 

To run a test case, execute ./Allclean && ./Allrun scripts  

NOTE: The test cases inside run/implicitContact might not run as they are under development.


# Post-processing using ParaView
    
The deformed configuration of the beam can be observed by 'warp by vector' using 'pointW' field.


# Note on the beamSolver:

1. The primary variables are displacement (W) and rotation (Theta) - 6 DOFs per cell

2. This solver does not use blockMesh. It has its own meshing utility. Can be found in applications/utilities/.

3. For starters, run the test cases in run/beamFoam/elastic/myCases/AllTestCasesInPaper/
 
   The test cases in this folder are also published in the paper:
   
A cell-centered finite volume formulation of geometrically exact Simo–Reissner beams with arbitrary initial curvatures, https://doi.org/10.1002/nme.6994
    
