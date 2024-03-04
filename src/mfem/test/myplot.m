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

fclose(fileID);

% ==========================================

M = size(X,2);

M = sqrt(M); 

X = reshape(X,M,M);
Y = reshape(Y,M,M);
vX = reshape(vX,M,M);
vY = reshape(vY,M,M);
vXX = reshape(vXX,M,M);
vYY = reshape(vYY,M,M);

figure
quiver(X,Y,vX,vY);
title("Approximation Velocity");

figure
quiver(X,Y,vXX,vYY);
title("Exact Velocity");




