function [] = myplot_melt(M, N, L, start)

% Input grid files
fileID = fopen('build/gridCellX.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridCellY.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

mid = floor(M/2) + 1;

v = VideoWriter('phase.avi','Motion JPEG AVI');
open(v);

fstruct1 = dir('build/*porosity*.dat');
fcell1 = struct2cell(fstruct1);

loops = numel(fstruct1)

h = figure;

set(gcf, 'Position',[100 100 1210 693])

for kk=1:loops-1

k = kk +start;

filename = strcat('build/porosity',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
porosity = fscanf(fileID, '%f', [1,Inf]);
porosity = reshape(porosity, M, N);

kn = k+1;

filename = strcat('build/porosity',string(kn));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
porosity2 = fscanf(fileID, '%f', [1,Inf]);
porosity2 = reshape(porosity2, M, N);


plot(porosity2(mid,:) - porosity(mid,:),pY(mid,:));
title(filename);
ylabel("Depth");
xlabel("Porosity");

pause
F= getframe(gcf);

writeVideo(v,F);

end
