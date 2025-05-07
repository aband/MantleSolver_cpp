function [] = myplot_velocity(M,N)

% Read grid files
fileID = fopen('build/gaussgridx.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gaussgridy.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

M = 3*M;

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

fstruct = dir('build/*darcyvelx*.dat');
fcell = struct2cell(fstruct);

loops = numel(fstruct);

for k = 1:loops

filename = strcat('build/darcyvelx',string(k));
filename = strcat(filename,'.dat')
fileID = fopen(filename, 'r');
vx = fscanf(fileID, '%f', [1,Inf]);
vx = reshape(vx, M, N);

filename = strcat('build/darcyvely',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, M, N);

subplot(2,2,1)
quiver(pX, pY, vx, vy);
title(["darcy velocity", num2str(k)])

subplot(2,2,3)
plot(vy(4,:), pY(4,:));

filename = strcat('build/stokesvelx',string(k));
filename = strcat(filename,'.dat')
fileID = fopen(filename, 'r');
vx = fscanf(fileID, '%f', [1,Inf]);
vx = reshape(vx, M, N);

filename = strcat('build/stokesvely',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, M, N);

subplot(2,2,2)
quiver(pX, pY, vx, vy);
title(["stokes velocity", num2str(k)])

subplot(2,2,4)
plot(vy(4,:), pY(4,:));

%set(gcf, 'Position',[50 50 1800 700]);
pause
G = getframe(gcf);

end
