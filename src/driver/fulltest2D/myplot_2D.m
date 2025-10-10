function [] = myplot_2D(M, N, start, folder)

% Input grid files
filename = strcat(folder, '/gridCellX.dat');
fileID = fopen(filename,'r');
pX = fscanf(fileID, '%f', [1,Inf]);

filename = strcat(folder, '/gridCellY.dat')
fileID = fopen(filename,'r');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

mid = floor(M/2) + 1;

v = VideoWriter('phase.avi','Motion JPEG AVI');
open(v);

filename = strcat(folder, '/*porosity*.dat');
fstruct1 = dir(filename);
fcell1 = struct2cell(fstruct1);

%fstruct2 = dir('build/*temperature*.dat');
%fcell2 = struct2cell(fstruct2);

loops = numel(fstruct1)

h = figure;

set(gcf, 'Position',[100 100 1210 693])

% Read grid files
filename = strcat(folder,'/gaussgridx.dat')
fileID = fopen(filename,'r');
ppX = fscanf(fileID, '%f', [1,Inf]);

filename = strcat(folder,'/gaussgridy.dat')
fileID = fopen(filename,'r');
ppY = fscanf(fileID, '%f', [1,Inf]);

MM = 3*M;

ppX = reshape(ppX, MM, N+1);
ppY = reshape(ppY, MM, N+1);


for kk=1:loops

k = kk + start;

filename = strcat(folder, '/porosity')
filename = strcat(filename,string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
porosity = fscanf(fileID, '%f', [1,Inf]);
porosity = reshape(porosity, M, N);


subplot(1,4,1)
surf(pX, pY, porosity);
view(2)
colormap turbo
shading interp
colorbar
%title("Porosity");
%ylabel("Depth");
title(["porosity"])

filename = strcat(folder, '/effvelx');
filename = strcat(filename,string(k));
filename = strcat(filename,'.dat')
fileID = fopen(filename, 'r');
vx = fscanf(fileID, '%f', [1,Inf]);
vx = reshape(vx, MM, N+1);

filename = strcat(folder, '/effvely');
filename = strcat(filename,string(k));
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

filename = strcat(folder, '/solidvelx');
filename = strcat(filename,string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vx = fscanf(fileID, '%f', [1,Inf]);
vx = reshape(vx, MM, N+1);

filename = strcat(folder, '/solidvely');
filename = strcat(filename,string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, MM, N+1);

subplot(1,4,3)
quiver(ppX, ppY, vx, vy);
%title(["Solid velocity", num2str(k)])
title(["Solid velocity"])

filename = strcat(folder, '/phase');
filename = strcat(filename, string(k));
filename = strcat(filename, '.dat');
fileID   = fopen(filename, 'r');
data     = fscanf(fileID, '%f', [1, Inf]);
data = reshape(data, M, N);

data = [data, ones(N,1);
        ones(1,M), 1];

subplot(1,4,4);
pcolor(data'-1)
colormap turbo
colorbar
clim([0,3.1])
title('phase')

time = 0.1*k*10;
mytitle= strcat('Time = ', string(time));
sgtitle(mytitle);

pause
G = getframe(gcf);

writeVideo(v,G);

end
