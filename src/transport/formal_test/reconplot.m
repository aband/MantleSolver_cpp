function [] = reconplot(M, N, mark, mytitle)

% Read grid files 

fileID = fopen('savedrun/gridreconx.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('savedrun/gridrecony.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, 3*M, 3*N);
pY = reshape(pY, 3*M, 3*N);

filename = strcat('savedrun/reconSol',string(mark));
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
