clc; clear

ux = h5read("build/velocity.cgns","/Base/Grid/Velocity/ux/ data");
uy = h5read("build/velocity.cgns","/Base/Grid/Velocity/uy/ data");

X  = h5read("build/velocity.cgns","/Base/Grid/GridCoordinates/CoordinateX/ data");
Y  = h5read("build/velocity.cgns","/Base/Grid/GridCoordinates/CoordinateY/ data");

X = X(1:end-1,1:end-1) + X(1:end-1, 2:end) + X(2:end,1:end-1) + X(2:end, 2:end);
X = X/4;

Y = Y(1:end-1,1:end-1) + Y(1:end-1, 2:end) + Y(2:end,1:end-1) + Y(2:end, 2:end);
Y = Y/4;

figure

quiver(X,Y,ux,uy);
title("Stokes Velocity");

figure
streamslice(X',Y',ux',uy',0.5);
title("Stokes Streamline");
