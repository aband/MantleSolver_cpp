function [] = reconplot(M, N, mark)

% Read grid files 

fileID = fopen('build/gridreconx.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridrecony.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, 3*M, 3*N);
pY = reshape(pY, 3*M, 3*N);

filename = strcat('build/reconSol',string(mark));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
sol = fscanf(fileID, '%f', [1,Inf]);
sol = reshape(sol, 3*M, 3*N);

surf(pX, pY, sol)
title(filename)
ylabel("y")
xlabel("x")
