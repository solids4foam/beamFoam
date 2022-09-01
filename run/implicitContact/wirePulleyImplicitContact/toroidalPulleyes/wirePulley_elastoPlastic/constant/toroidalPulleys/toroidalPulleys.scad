$fn = 100;

translate([0, 0, 0])
rotate([0, 0, 0])
difference() {cylinder(10, 18.0013, 18.0013, true);
rotate_extrude(convexity=10) translate([10018, 0, 0]) circle(10000);}

translate([60, 0, 0])
rotate([0, 0, 0])
difference() {cylinder(10, 19.5013, 19.5013, true);
rotate_extrude(convexity=10) translate([10019.5, 0, 0]) circle(10000);}

translate([125, -10, 0])
rotate([0, 0, 0])
difference() {cylinder(10, 25.5013, 25.5013, true);
rotate_extrude(convexity=10) translate([10025.5, 0, 0]) circle(10000);}

