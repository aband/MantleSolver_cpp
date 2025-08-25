function [] = myplot_cut(M,N, cut)

% Print porosity and velocity at the given time stamp

fileID = fopen('build/gridCellX.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridCellY.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

filename = strcat('build/porosity', string(cut));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
poro = fscanf(fileID, '%f', [1,Inf]);
poro = reshape(poro, M, N);

subplot(1,2,2)
plot(poro(2,:), pY(2,:));
title(filename)
ylabel("Depth");
axis([])
%axis equal
