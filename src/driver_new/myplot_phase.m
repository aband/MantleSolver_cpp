function [] = phaseplot(M, N, mark, folder, field)

%Read grid files

filename = strcat(folder, '/build/recongridx');
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
gridx = fscanf(fileID, '%f', [1,Inf]);
gridx = reshape(gridx, 3*M, 3*N);

filename = strcat(folder, '/build/recongridy');
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
gridy = fscanf(fileID, '%f', [1,Inf]);
gridy = reshape(gridy, 3*M, 3*N);

filename = strcat(folder, '/build/');
filename = strcat(filename, field);
filename = strcat(filename,string(mark));
filename = strcat(filename,'.dat')
fileID = fopen(filename, 'r');
sol = fscanf(fileID, '%f', [1,Inf]);
sol = reshape(sol, 3*M, 3*N);


figure
s = surf(gridx, gridy, sol)
s.EdgeColor = 'none';
colormap(turbo)
