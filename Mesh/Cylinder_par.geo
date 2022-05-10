// Gmsh project created on Fri Feb 18 13:17:48 2022
d=DefineNumber[0.1, Name "Parameters/d"];
r=DefineNumber[d/2, Name "Parameters/r"];
cos=DefineNumber[0.70710678*d*0.5, Name "Parameters/cos"];
r_inn=DefineNumber[0.07, Name "Parameters/r_inn"];
H=DefineNumber[0.41, Name "Parameters/H"];
L_front=DefineNumber[0.2, Name "Parameters/L_front"];
L_back=DefineNumber[2, Name "Parameters/L_back"];
SetFactory("OpenCASCADE");
//+
Point(1) = {-L_front, -H, 0, 1.0};
//+
Point(2) = {L_back, -H, 0, 1.0};
//+
Point(3) = {L_back, H, 0, 1.0};
//+
Point(4) = {-L_front, H, 0, 1.0};
//+
Point(5) = {-cos-r_inn, -H, 0, 1.0};
//+
Point(6) = {cos+r_inn, -H, 0, 1.0};
//+
Point(7) = {-cos-r_inn, H, 0, 1.0};
//+
Point(8) = {cos+r_inn, H, 0, 1.0};
//+
Point(9) = {L_back, cos+r_inn, 0, 1.0};
//+
Point(10) = {L_back, -cos-r_inn, 0, 1.0};
//+
Point(11) = {-L_front, cos+r_inn, 0, 1.0};
//+
Point(12) = {-L_front, -cos-r_inn, 0, 1.0};
//+
Point(13) = {cos+r_inn, cos+r_inn, 0, 1.0};
//+
Point(14) = {cos+r_inn, -cos-r_inn, 0, 1.0};
//+
Point(15) = {-cos-r_inn, cos+r_inn, 0, 1.0};
//+
Point(16) = {-cos-r_inn, -cos-r_inn, 0, 1.0};
//+
Point(17) = {0, 0, 0, 1.0};
//+
Point(18) = {cos, cos, 0, 1.0};
//+
Point(19) = {cos, -cos, 0, 1.0};
//+
Point(20) = {-cos, cos, 0, 1.0};
//+
Point(21) = {-cos, -cos, 0, 1.0};
//+
Circle(1) = {13, 17, 15};
//+
Circle(2) = {15, 17, 16};
//+
Circle(3) = {16, 17, 14};
//+
Circle(4) = {14, 17, 13};
//+
Line(5) = {1, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {2, 6};
//+
Line(8) = {2, 10};
//+
Line(9) = {10, 9};
//+
Line(10) = {3, 9};
//+
Line(11) = {3, 8};
//+
Line(12) = {8, 7};
//+
Line(13) = {4, 7};
//+
Line(14) = {4, 11};
//+
Line(15) = {11, 12};
//+
Line(16) = {1, 12};
//+
Line(17) = {12, 16};
//+
Line(18) = {11, 15};
//+
Line(19) = {9, 13};
//+
Line(20) = {10, 14};
//+
Line(21) = {7, 15};
//+
Line(22) = {8, 13};
//+
Line(23) = {6, 14};
//+
Line(24) = {5, 16};
//+
Circle(25) = {18, 17, 20};
//+
Circle(26) = {20, 17, 21};
//+
Circle(27) = {21, 17, 19};
//+
Circle(28) = {19, 17, 18};
//+
Line(29) = {18, 13};
//+
Line(30) = {19, 14};
//+
Line(31) = {21, 16};
//+
Line(32) = {20, 15};
//+
Curve Loop(1) = {14, 18, -21, -13};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {21, -1, -22, 12};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {22, -19, -10, 11};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {4, -19, -9, 20};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {20, -23, -7, 8};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {3, -23, -6, 24};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {17, -24, -5, 16};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {17, -2, -18, 15};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {32, -1, -29, 25};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {29, -4, -30, 28};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {27, 30, -3, -31};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {26, 31, -2, -32};
//+
Plane Surface(12) = {12};
//+
Physical Curve("Inlet", 25) = {14, 15, 16};
//+
Physical Curve("Outlet", 26) = {10, 9, 8};
//+
Physical Curve("Limits", 27) = {13, 12, 11, 5, 6, 7};
//+
Physical Curve("Outer_Cylinder", 28) = {2, 1, 4, 3};
//+
Physical Curve("Cylinder", 33) = {26, 25, 28, 27};
//+
Physical Surface("Fluid", 29) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
//+
Recombine Surface {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
//+
Transfinite Curve {15, 2, 1, 4, 3, 12, 9, 6} = 21 Using Progression 1;
//+
Transfinite Curve {14, 21, 22, 10, 8, 23, 24, 16} = 21 Using Progression 0.95;
//+
Transfinite Curve {13, 18, 17, 5} = 11 Using Progression 0.95;
//+
Transfinite Curve {11, 19, 20, 7} = 71 Using Progression 0.98;
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Transfinite Surface {5};
//+
Transfinite Surface {6};
//+
Transfinite Surface {7};
//+
Transfinite Surface {8};
//+
Transfinite Surface {9};
//+
Transfinite Surface {10};
//+
Transfinite Surface {11};
//+
Transfinite Surface {12};
//+
Transfinite Curve {2, 1, 4, 3, 27, 28, 25, 26} = 21 Using Progression 1;
//+
Transfinite Curve {29, 30, 31, 32} = 40 Using Progression 1.1;



//+
Physical Point("P_inlet", 34) = {11, 12};
//+
Physical Point("P_outlet", 35) = {9, 10};
//+
Physical Point("P_limits", 36) = {4, 7, 8, 3, 2, 6, 5, 1};
//+
Physical Point("P_cylinder", 37) = {20, 18, 19, 21};
