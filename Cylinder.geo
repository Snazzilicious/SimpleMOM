

rad = 1.0;
Len = 10.0;
scale = 0.1;

Point(1) = {0.0, 0.0, -Len/2.0, scale};

Point(2) = {rad, 0.0, -Len/2.0, scale};
Point(3) = {0.0, rad, -Len/2.0, scale};
Point(4) = {-rad, 0.0, -Len/2.0, scale};
Point(5) = {0.0, -rad, -Len/2.0, scale};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};


Point(6) = {0.0, 0.0, Len/2.0, scale};

Point(7) = {rad, 0.0, Len/2.0, scale};
Point(8) = {0.0, rad, Len/2.0, scale};
Point(9) = {-rad, 0.0, Len/2.0, scale};
Point(10) = {0.0, -rad, Len/2.0, scale};

Circle(5) = {7,6,8};
Circle(6) = {8,6,9};
Circle(7) = {9,6,10};
Circle(8) = {10,6,7};


Line(9) = {2,7};
Line(10) = {3,8};
Line(11) = {4,9};
Line(12) = {5,10};


Curve Loop(1) = {1,10,-5,-9};
Curve Loop(2) = {2,11,-6,-10};
Curve Loop(3) = {3,12,-7,-11};
Curve Loop(4) = {4,9,-8,-12};


For i In {1:4}
  Surface(i) = {i};
EndFor
