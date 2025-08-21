function [] = reconplot(M, N, mark)

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


filename = strcat('build/reconSol',string(mark));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
sol2 = fscanf(fileID, '%f', [1,Inf]);
sol2 = reshape(sol2, 3*M, 3*N);

figure
s = surf(pX, pY, sol)
s.EdgeColor = 'none';
title('exact')
ylabel("y")
xlabel("x")
colormap(turbo)
colorbar
%caxis([0,1])

figure
s = surf(pX, pY, sol2)
s.EdgeColor = 'none';
title('reconstructed')
ylabel("y")
xlabel("x")
colormap(turbo)
colorbar
%caxis([0,1])
