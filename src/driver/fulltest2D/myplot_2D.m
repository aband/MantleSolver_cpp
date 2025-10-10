function [] = myplot_2D(M, N, start)

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

fstruct2 = dir('build/*temperature*.dat');
fcell2 = struct2cell(fstruct2);

loops = numel(fstruct1)

h = figure;

set(gcf, 'Position',[100 100 1210 693])

% Read grid files
fileID = fopen('build/gaussgridx.dat','r');
ppX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gaussgridy.dat','r');
ppY = fscanf(fileID, '%f', [1,Inf]);

MM = 3*M;

ppX = reshape(ppX, MM, N+1);
ppY = reshape(ppY, MM, N+1);


for kk=1:loops

k = kk + start;

filename = strcat('build/porosity',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
porosity = fscanf(fileID, '%f', [1,Inf]);
porosity = reshape(porosity, M, N);

subplot(1,4,1)
surf(pX, pY, porosity);
colormap turbo
shading interp
%title("Porosity");
ylabel("Depth");
xlabel("Porosity");

filename = strcat('build/effvelx',string(k));
filename = strcat(filename,'.dat')
fileID = fopen(filename, 'r');
vx = fscanf(fileID, '%f', [1,Inf]);
vx = reshape(vx, MM, N+1);

filename = strcat('build/effvely',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, MM, N+1);

subplot(1,4,2)
quiver(ppX, ppY, vx, vy);
%title(["Effective velocity", num2str(k)])
title(["Effective velocity"])

%{
filename = strcat('build/phasevelx',string(k));
filename = strcat(filename,'.dat')
fileID = fopen(filename, 'r');
vx = fscanf(fileID, '%f', [1,Inf]);
vx = reshape(vx, MM, N+1);

filename = strcat('build/phasevely',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, MM, N+1);

subplot(1,3,3)
quiver(ppX, ppY, vx, vy);
%title(["Phase averaged velocity", num2str(k)])
title(["Phase averaged velocity"])

%subplot(2,3,5)
%plot(vy(4,:), pY(4,:), 'LineWidth', 3);
%xlim([-1,1])

%}

filename = strcat('build/solidvelx',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vx = fscanf(fileID, '%f', [1,Inf]);
vx = reshape(vx, MM, N+1);

filename = strcat('build/solidvely',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, MM, N+1);

subplot(1,4,3)
quiver(ppX, ppY, vx, vy);
%title(["Solid velocity", num2str(k)])
title(["Solid velocity"])

filename = strcat('build/phase', string(k));
filename = strcat(filename, '.dat');
fileID   = fopen(filename, 'r');
data     = fscanf(fileID, '%f', [1, Inf]);
data = reshape(data, M, N);

subplot(1,6,5);
surf(pX, pY, phase)

pause

end


