function [] = shape(N)

fileID = fopen('build/referenceX.dat','r');
X = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/referenceY.dat','r');
Y = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/valx.dat','r');
vX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/valy.dat','r');
vY = fscanf(fileID, '%f', [1,Inf]);

X = reshape(X,N,N);
Y = reshape(Y,N,N);
vX = reshape(vX,N,N);
vY = reshape(vY,N,N);

figure
quiver(X,Y,vX,vY);
title("Shape Function");
