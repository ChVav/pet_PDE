lc = 0.05;

Point(1) = {0,0,0,lc};
Point(2) = {0,1,0,lc};
Point(3) = {1,1,0,lc};
Point(4) = {1,0,0,lc};
//+
Line(1) = {2, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 1};
//+
Line(4) = {1, 2};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Curve("boundary_up", 5) = {1};
//+
Physical Curve("boundary_low", 6) = {3};
//+
Physical Curve("boundary_left", 7) = {4};
//+
Physical Curve("boundary_right", 8) = {2};
