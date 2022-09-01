# Details of this test case: 

This perpendicular beams test case is setup to test the ability of pointContact 
when the contact occurs exactly at the face of the cell.

Observations as of today (31st August 2022):

When the contact is exactly at thhe face of the cell, the solver does not converge 
per time step and the contact point in the search algorithm jumps from one cell to the adjacent
cell after every iteration in a time step (thus not converging)

An attempt was made to forcefully make the 
beam weights 0.5,  if zeta (location of point contact) > 0.998 || < -0.998
but this does not solve the problem.

Follow up with ZT and work on this test case.

