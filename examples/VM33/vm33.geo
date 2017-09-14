//+
lc = DefineNumber[ 0.3, Name "Parameters/lc" ];
//+
Point(1) = {0, 0, 0, lc};
//+
Point(2) = {1, 0, 0, lc};
//+
Point(3) = {3, 0, 0, lc};
//+
Point(4) = {0, 1, 0, lc};
//+
Point(5) = {0, 3, 0, lc};
//+
Circle(1) = {2, 1, 4};
//+
Circle(2) = {3, 1, 5};
//+
Line(3) = {2, 3};
//+
Line(4) = {4, 5};
//+
Line Loop(1) = {1, 4, -2, -3};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, 0.4} {
  Surface{1}; 
}
//+
Physical Surface("outter") = {21};
//+
Physical Surface("zeroFux") = {26, 1, 13};
//+
Physical Surface("zeroFux") += {25, 17};
//+
Physical Volume("mat1") = {1};
