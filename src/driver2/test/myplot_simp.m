function [] = myplot_1d(M, N, name)

% Read grid files 

fileID = fopen('build/gridCellX.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridCellY.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

%v = VideoWriter('video.avi','Motion JPEG AVI');
%open(v);

fullname = strcat('build/', name);
fullname = strcat(fullname, '*.dat');
fstruct1 = dir(fullname);
fcell1 = struct2cell(fstruct1);

loops = numel(fstruct1)

h = figure;

for k=1:loops

% Porosity
filename = strcat('build/', name);
filename = strcat(filename, string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
sol = fscanf(fileID, '%f', [1,Inf]);
sol = reshape(sol, M, N);

%surf(pX, pY, sol)

plot(sol(2,:), pY(2,:));
title(name)
ylabel("Depth");

set(gcf, 'Position',[50 50 1800 700]);

%pause
G = getframe(gcf);

%writeVideo(v,G);
end

%close(v);
