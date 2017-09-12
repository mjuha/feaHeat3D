//+
SetFactory("OpenCASCADE");
Cylinder(1) = {0, 0, 0, 0, 0.01, 0, 0.005, 2*Pi};
//+
Physical Volume("mat1") = {1};
//+
Physical Surface("convLoad") = {2, 1, 3};
