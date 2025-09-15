function [] = myplot_1d(M, N, name, folder)

% Read grid files 

fileID = fopen('build/gridCellX.dat','r');
%fileID = fopen('record/gridCellX.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridCellY.dat','r');
%fileID = fopen('record/gridCellY.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, M, N)%;
pY = reshape(pY, M, N);

v = VideoWriter('video.avi','Motion JPEG AVI');
open(v);

fullname = strcat(folder, name);
fullname = strcat(fullname, '*.dat');
fstruct1 = dir(fullname);
fcell1 = struct2cell(fstruct1);

loops = numel(fstruct1)

h = figure;

set(gcf, 'Position',[50 50 1800 700]);

for k=1:loops

% Porosity
filename = strcat(folder, name);
filename = strcat(filename, string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
sol = fscanf(fileID, '%f', [1,Inf]);
sol = reshape(sol, M, N);

subplot(1,2,1)
surf(pX, pY, sol)

subplot(1,2,2)
plot(sol(2,:), pY(2,:));
title(filename)
ylabel("Depth");

pause
G = getframe(gcf);

writeVideo(v,G);
end

close(v);
