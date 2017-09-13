//+
SetFactory("OpenCASCADE");
Cylinder(1) = {0, 0, 0, 0, 0.1, 0, 0.025, 2*Pi};
//+
Physical Volume("mat1") = {1};
//+
Physical Surface("convLoad") = {2, 1, 3};
//+
Mesh.CharacteristicLengthFactor = 0.9;
