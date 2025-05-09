function [] = myplot_velocity(M,N)

% Read grid files
fileID = fopen('build/gaussgridx.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gaussgridy.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

M = 3*M;

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

fstruct = dir('build/*effvel_pause*.dat');
fcell = struct2cell(fstruct);

loops = numel(fstruct);

for k = 1:loops

filename = strcat('build/effvel_pause',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, M, N);

subplot(2,3,1)
quiver(pX, pY, 0.0*vy, vy);
title(["Effective velocity", num2str(k)])

subplot(2,3,4)
plot(vy(4,:), pY(4,:));

filename = strcat('build/phasevel_pause',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, M, N);

subplot(2,3,2)
quiver(pX, pY, 0.0*vy, vy);
title(["Phase averaged velocity", num2str(k)])

subplot(2,3,5)
plot(vy(4,:), pY(4,:));

filename = strcat('build/solidvel_pause',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, M, N);

subplot(2,3,3)
quiver(pX, pY, 0.0*vy, vy);
title(["Solid velocity", num2str(k)])

subplot(2,3,6)
plot(vy(4,:), pY(4,:));

%set(gcf, 'Position',[50 50 1800 700]);
pause
G = getframe(gcf);

end
