# Test case details:

An initially curved arc is subjected to a point load at the crown (apex).

Details can be found in the paper: A cell-centered finite volume formulation of geometrically exact Simo–Reissner beams with arbitrary initial curvatures, https://doi.org/10.1002/nme.6994

Reference configuration - The beam is straight. 

# Initial configuration

makeCircularArc - sets the initial configuration of the beam from straight to an arc and updates the reference displacement (W) and rotation (Theta) fields.

You can visualise the deformed configuration in paraview by warping the vector using 'pointW' field

# Location of the point Force

The location of the force is set in constant/pointForces 
E.g. If the number os beam Segments is 41 (odd number so that crown location coincides with the cell centre of 20th cell)
Hence, the location where the pointForce is applied is: ( 20 0) where 0 indicates the cell centre.



  
