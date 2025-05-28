function [] = myplot_velocity(M,N)

fileID = fopen('build/gridCellX.dat','r');
cX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridCellY.dat','r');
cY = fscanf(fileID, '%f', [1,Inf]);

cX = reshape(cX, M, N);
cY = reshape(cY, M, N);

% Read grid files
fileID = fopen('build/gaussgridx.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gaussgridy.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

M = 3*M;
N = N+1;

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

fstruct = dir('build/*darcyvelx*.dat');
fcell = struct2cell(fstruct);

loops = numel(fstruct);

%set(gcf, 'Position',[50 50 1800 700]);
for k = 1:loops

figure

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
title(["darcy velocity"])

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
title(["stokes velocity"])

subplot(2,2,4)
plot(vy(4,:), pY(4,:));

figure
filename = strcat('build/qs',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
qs = fscanf(fileID, '%f', [1,Inf]);
qs = reshape(qs, 2, 100);

subplot(1,2,1)
plot(qs(1,:), cY(1,:));

filename = strcat('build/qf',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
qf = fscanf(fileID, '%f', [1,Inf]);
qf = reshape(qf, 2, 100);

subplot(1,2,2)
plot(qf(1,:), cY(1,:));

end
