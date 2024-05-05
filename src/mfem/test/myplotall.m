function [] = myplotall(M,N)

fileID = fopen('build/gridX.dat','r');
X = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridY.dat','r');
Y = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/vx.dat','r');
vX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/vy.dat','r');
vY = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/ux.dat','r');
uX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/uy.dat','r');
uY = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/porosity.dat','r');
poro = fscanf(fileID, '%f', [1,Inf]);

%{
fileID = fopen('sample/gridX.dat','r');
X = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('sample/gridY.dat','r');
Y = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('sample/vx.dat','r');
vX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('sample/vy.dat','r');
vY = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('sample/ux.dat','r');
uX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('sample/uy.dat','r');
uY = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('sample/porosity.dat','r');
poro = fscanf(fileID, '%f', [1,Inf]);
%}

fclose(fileID);

% =========================================

X = reshape(X,M,N);
Y = reshape(Y,M,N);
vX = reshape(vX,M,N);
vY = reshape(vY,M,N);
uX = reshape(uX,M,N);
uY = reshape(uY,M,N);
poro = reshape(poro,M,N);

figure
quiver(X,Y,vX,vY);
title("Stokes Velocity");

figure
streamslice(X',Y',vX',vY',0.5);
title("Stokes Streamline");

figure
quiver(X,Y,uX,uY);
title("Darcy Velocity");

figure
streamslice(X',Y',uX',uY',0.5);
title("Darcy Streamline");

figure
contour(X,Y,poro,20);
title("Porosity");
