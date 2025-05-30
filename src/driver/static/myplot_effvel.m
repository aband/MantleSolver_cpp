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

filename = ('build/exactv.dat');
fileID = fopen(filename,'r')
exact = fscanf(fileID, '%f', [1,Inf]);
exact = reshape(exact, M, N);

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
hold on
plot(vy(4,:), pY(4,:));
plot(-1*exact(4,:), pY(4,:),'+');
hold off

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
hold on
plot(vy(4,:), pY(4,:));
plot(exact(4,:), pY(4,:), '+');
hold off

figure
filename = strcat('build/qs',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
qs = fscanf(fileID, '%f', [1,Inf]);
qs = reshape(qs, M/3, N-1);

subplot(1,4,1)
plot(qs(1,:), cY(1,:));
title("qs")

filename = strcat('build/qf',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
qf = fscanf(fileID, '%f', [1,Inf]);
qf = reshape(qf, M/3, N-1);

subplot(1,4,2)
plot(qf(1,:), cY(1,:));
title("qf")

filename = strcat('build/rawstokesq',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
rqs = fscanf(fileID, '%f', [1,Inf]);
rqs = reshape(rqs, M/3, N-1);

subplot(1,4,3)
plot(rqs(1,:), cY(1,:));
title("raw data qs")

filename = strcat('build/rawdarcyq',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
rqf = fscanf(fileID, '%f', [1,Inf]);
rqf = reshape(rqf, M/3, N-1);

subplot(1,4,4)
plot(rqf(1,:), cY(1,:));
title("raw data qf")

end
