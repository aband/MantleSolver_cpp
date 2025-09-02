function [] = myplot_velocity(M,N, mark, folder)

fileID = strcat(folder, '/build/gaussgridx.dat')
fileID = fopen(fileID,'r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = strcat(folder, '/build/gaussgridy.dat');
fileID = fopen(fileID,'r');
pY = fscanf(fileID, '%f', [1,Inf]);

M = 3*M;

pY
size(pX)
size(pY)

pX = reshape(pX, M, N+1);
pY = reshape(pY, M, N+1);

% ====================================================
filename = strcat(folder, '/build/effvelx');
filename = strcat(filename, string(mark));
filename = strcat(filename,'.dat')
fileID = fopen(filename, 'r');
vx = fscanf(fileID, '%f', [1,Inf]);
vx = reshape(vx, M, N+1);

filename = strcat(folder, '/build/effvely');
filename = strcat(filename, string(mark));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, M, N+1);

subplot(2,3,1)
quiver(pX, pY, vx, vy);
title(["Effective velocity", num2str(mark)])

subplot(2,3,4)
plot(vy(4,:), pY(4,:));

% ====================================================
filename = strcat(folder, '/build/phasevelx');
filename = strcat(filename, string(mark));
filename = strcat(filename,'.dat')
fileID = fopen(filename, 'r');
vx = fscanf(fileID, '%f', [1,Inf]);
vx = reshape(vx, M, N+1);

filename = strcat(folder, '/build/phasevely');
filename = strcat(filename, string(mark));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, M, N+1);

subplot(2,3,2)
quiver(pX, pY, vx, vy);
title(["Phase averaged velocity", num2str(mark)])

subplot(2,3,5)
plot(vy(4,:), pY(4,:));
xlim([-1,1])

% ====================================================
filename = strcat(folder,'/build/solidvely');
filename = strcat(filename,string(mark));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, M, N+1);

subplot(2,3,3)
quiver(pX, pY, vx, vy);
title(["Solid velocity", num2str(mark)])

subplot(2,3,6)
plot(vy(4,:), pY(4,:));
