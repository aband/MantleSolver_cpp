function [] = myplot_1d(M, N)

% Read grid files 

fileID = fopen('build/gridCellX.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridCellY.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

v = VideoWriter('test.avi','Motion JPEG AVI');
open(v);

fstruct1 = dir('build/*sol*.dat');
fcell1 = struct2cell(fstruct1);

loops = numel(fstruct1)

h = figure;

set(gcf, 'Position',[50 50 1800 700]);

for k=1:loops

% Porosity
filename = strcat('build/sol',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
sol = fscanf(fileID, '%f', [1,Inf]);
sol = reshape(sol, M, N);

surf(pX, pY, sol)
title(filename)
ylabel("Depth");

pause
G = getframe(gcf);

writeVideo(v,G);
end

close(v);
