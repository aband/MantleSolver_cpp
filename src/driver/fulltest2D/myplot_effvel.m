function [] = myplot_velocity(M,N)

% Read grid files
fileID = fopen('build/gaussgridx.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gaussgridy.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

M = 3*M;

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

fstruct = dir('build/*effvelx*.dat');
fcell = struct2cell(fstruct);

loops = numel(fstruct);

v = VideoWriter('edgevel.avi','Motion JPEG AVI');
open(v);

h = figure;

set(gcf, 'Position',[100 100 1210 693])

for k = 1:loops

filename = strcat('build/effvelx',string(k));
filename = strcat(filename,'.dat')
fileID = fopen(filename, 'r');
vx = fscanf(fileID, '%f', [1,Inf]);
vx = reshape(vx, M, N);

filename = strcat('build/effvely',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, M, N);

subplot(1,3,1)
quiver(pX, pY, vx, vy);
%title(["Effective velocity", num2str(k)])
title(["Effective velocity"])

%subplot(2,3,4)
%plot(vy(4,:), pY(4,:), 'LineWidth', 3);

filename = strcat('build/phasevelx',string(k));
filename = strcat(filename,'.dat')
fileID = fopen(filename, 'r');
vx = fscanf(fileID, '%f', [1,Inf]);
vx = reshape(vx, M, N);

filename = strcat('build/phasevely',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, M, N);

subplot(1,3,2)
quiver(pX, pY, vx, vy);
%title(["Phase averaged velocity", num2str(k)])
title(["Phase averaged velocity"])

%subplot(2,3,5)
%plot(vy(4,:), pY(4,:), 'LineWidth', 3);
%xlim([-1,1])

filename = strcat('build/solidvely',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, M, N);

subplot(1,3,3)
quiver(pX, pY, vx, vy);
%title(["Solid velocity", num2str(k)])
title(["Solid velocity"])

%subplot(2,3,6)
%plot(vy(4,:), pY(4,:), 'LineWidth', 3);

time = 0.1*k*10;

mytitle = strcat('Time = ', string(time));

sgtitle(mytitle);

%set(gcf, 'Position',[50 50 1800 700]);
pause
G = getframe(gcf);

writeVideo(v,G);

end
