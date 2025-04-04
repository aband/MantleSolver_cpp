function [] = myplot_velocity(M,N)

% Read grid files
fileID = fopen('build/gaussgridx.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gaussgridy.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

M = 3*M;

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

fstruct = dir('build/*effvelxHD*.dat');
fcell = struct2cell(fstruct);

loops = numel(fstruct);

for k = 1:loops

filename = strcat('build/effvelxHD',string(k));
filename = strcat(filename,'.dat')
fileID = fopen(filename, 'r');
vxhd = fscanf(fileID, '%f', [1,Inf]);
vxhd = reshape(vxhd, M, N);

filename = strcat('build/effvelyHD',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vyhd = fscanf(fileID, '%f', [1,Inf]);
vyhd = reshape(vyhd, M, N);

filename = strcat('build/effvelxCD',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vxcd = fscanf(fileID, '%f', [1,Inf]);
vxcd = reshape(vxcd, M, N);

filename = strcat('build/effvelyCD',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vycd = fscanf(fileID, '%f', [1,Inf]);
vycd = reshape(vycd, M, N);

subplot(1,2,1)
quiver(pX, pY, vxhd, vyhd);
title(["Effective velocity HD", num2str(k)])

subplot(1,2,2)
quiver(pX, pY, vxcd, vycd);
title(["Effective velocity CD", num2str(k)])

set(gcf, 'Position',[50 50 1800 700]);

G = getframe(gcf);

end
