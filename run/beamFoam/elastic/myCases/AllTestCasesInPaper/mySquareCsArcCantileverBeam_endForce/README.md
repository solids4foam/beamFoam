# Test case details:

An initially 45 deg curved cantilever is subjected to a point load at the right end.

Details can be found in the paper: A cell-centered finite volume formulation of geometrically exact Simo–Reissner beams with arbitrary initial curvatures, https://doi.org/10.1002/nme.6994

Reference configuration - The beam is straight. 

# Initial configuration

makeArc - sets the initial configuration of the beam from straight to an 45 deg arc and updates the reference displacement (W) and rotation (Theta) fields.

You can visualise the deformed configuration in paraview by warping the vector using 'pointW' field




  
