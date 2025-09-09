function [] = myplot_velocity(M,N,folder1,folder2)

% Read grid files
%fileID = fopen('build/gaussgridx.dat','r');
filename = strcat(folder2, '/gaussgridx');
filename = strcat('/', filename);
filename = strcat(folder1, filename);
filename = strcat(filename, '.dat');
fileID = fopen(filename,'r');
pX = fscanf(fileID, '%f', [1,Inf]);

filename = strcat(folder2, '/gaussgridy');
filename = strcat('/', filename);
filename = strcat(folder1, filename);
filename = strcat(filename, '.dat');
fileID = fopen(filename,'r');
pY = fscanf(fileID, '%f', [1,Inf]);

M = 3*M;

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

filename = strcat(folder2, '/*effvelx*.dat');
filename = strcat('/', filename);
filename = strcat(folder1, filename);
fstruct = dir(filename);
fcell = struct2cell(fstruct);

loops = numel(fstruct);

v = VideoWriter('edgevel.avi','Motion JPEG AVI');
open(v);

h = figure;

set(gcf, 'Position',[100 100 1210 693])

for k = 1:loops

filename = strcat(folder2, '/effvelx');
filename = strcat(filename,string(k));
filename = strcat(filename,'.dat');
filename = strcat('/',filename);
filename = strcat(folder1, filename);
fileID = fopen(filename, 'r');
vx = fscanf(fileID, '%f', [1,Inf]);
vx = reshape(vx, M, N);


filename = strcat(folder2, '/effvely');
filename = strcat(filename,string(k));
filename = strcat(filename,'.dat');
filename = strcat('/',filename);
filename = strcat(folder1, filename);
fileID = fopen(filename, 'r');
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, M, N);

subplot(2,3,1)
quiver(pX, pY, vx, vy);
title(["Effective velocity", num2str(k)])

subplot(2,3,4)
plot(vy(4,:), pY(4,:));

filename = strcat(folder2, '/phasevelx');
filename = strcat(filename,string(k));
filename = strcat(filename,'.dat');
filename = strcat('/',filename);
filename = strcat(folder1, filename);
fileID = fopen(filename, 'r');
vx = fscanf(fileID, '%f', [1,Inf]);
vx = reshape(vx, M, N);

filename = strcat(folder2, '/phasevely');
filename = strcat(filename,string(k));
filename = strcat(filename,'.dat');
filename = strcat('/',filename);
filename = strcat(folder1, filename);
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, M, N)

subplot(2,3,2)
quiver(pX, pY, vx, vy);
title(["Phase averaged velocity", num2str(k)])

subplot(2,3,5)
plot(vy(4,:), pY(4,:));
xlim([-1e-4,1e-4])

filename = strcat(folder2, '/solidvely');
filename = strcat(filename,string(k));
filename = strcat(filename,'.dat');
filename = strcat('/',filename);
filename = strcat(folder1, filename);
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, M, N);

subplot(2,3,3)
quiver(pX, pY, vx, vy);
title(["Solid velocity", num2str(k)])

subplot(2,3,6)
plot(vy(4,:), pY(4,:));

%set(gcf, 'Position',[50 50 1800 700]);
pause
G = getframe(gcf);

writeVideo(v,G);

end
