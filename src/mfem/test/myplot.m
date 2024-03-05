clc; clear

% Read mesh points
% Read form gridX and gridY text file

fileID = fopen('build/gridX.dat','r');
X = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridY.dat','r');
Y = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/aprxVx.dat','r');
vX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/aprxVy.dat','r');
vY = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/exctVx.dat','r');
vXX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/exctVy.dat','r');
vYY = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/porosity.dat','r');
poro = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

% ==========================================

%M = size(X,2)
M = 80;
N = 40;
%M = sqrt(M); 

X = reshape(X,M,N);
Y = reshape(Y,M,N);
vX = reshape(vX,M,N);
vY = reshape(vY,M,N);
vXX = reshape(vXX,M,N);
vYY = reshape(vYY,M,N);
poro = reshape(poro,M,N);

figure
quiver(X,Y,vX,vY);
title("Approximation Velocity");

figure
quiver(X,Y,vXX,vYY);
title("Exact Velocity");

figure
contour(X,Y,poro);
title("Porosity");


