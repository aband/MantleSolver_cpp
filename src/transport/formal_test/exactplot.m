function [] = reconplot(M, N, mark, mytitle)

% Read grid files

fileID = fopen('build/exactgridx.dat');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/exactgridy.dat');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, 3*M, 3*N);
pY = reshape(pY, 3*M, 3*N);

filename = strcat('build/exactsol',string(mark));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
sol = fscanf(fileID, '%f', [1,Inf]);
sol = reshape(sol, 3*M, 3*N);

s = surf(pX, pY, sol)
s.EdgeColor = 'none';
title(mytitle)
ylabel("y")
xlabel("x")
colormap(turbo)
colorbar
caxis([0,1])
