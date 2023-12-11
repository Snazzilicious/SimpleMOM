

rad = 1.0;
scale = 0.1;

Point(1) = {0.0, 0.0, 0.0, scale};

Point(2) = {rad, 0.0, 0.0, scale};
Point(3) = {0.0, rad, 0.0, scale};
Point(4) = {-rad, 0.0, 0.0, scale};
Point(5) = {0.0, -rad, 0.0, scale};

Point(6) = {0.0, 0.0, rad, scale};
Point(7) = {0.0, 0.0, -rad, scale};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Circle(5) = {6,1,3};
Circle(6) = {3,1,7};
Circle(7) = {7,1,5};
Circle(8) = {5,1,6};

Circle(9) = {2,1,7};
Circle(10) = {7,1,4};
Circle(11) = {4,1,6};
Circle(12) = {6,1,2};


Curve Loop(1) = {1,-5,12};
Curve Loop(2) = {1,6,-9};
Curve Loop(3) = {2,-10,-6};
Curve Loop(4) = {2,11,5};
Curve Loop(5) = {3,8,-11};
Curve Loop(6) = {3,-7,10};
Curve Loop(7) = {4,9,7};
Curve Loop(8) = {4,-12,-8};


For i In {1:8}
  Surface(i) = {i} In Sphere{1};
EndFor
