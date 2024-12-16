clc; clear

ux = h5read("build/stokes.cgns","/Base/Grid/Velocity/ux/ data");
uy = h5read("build/stokes.cgns","/Base/Grid/Velocity/uy/ data");

vx = h5read("build/darcy.cgns","/Base/Grid/Velocity/ux/ data");
vy = h5read("build/darcy.cgns","/Base/Grid/Velocity/uy/ data");

X  = h5read("build/stokes.cgns","/Base/Grid/GridCoordinates/CoordinateX/ data");
Y  = h5read("build/stokes.cgns","/Base/Grid/GridCoordinates/CoordinateY/ data");

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

figure
quiver(X,Y,vx,vy);
title("Darcy Velocity");

figure
streamslice(X',Y',(abs(vx')>10e-6).*vx',(abs(vy')>10e-6).*vy',0.5);
title("Darcy Streamline");
