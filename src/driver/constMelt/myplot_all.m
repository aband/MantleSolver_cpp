function [] = myplot_all(M, N, folder)

filename = strcat(folder, '/gridCellX.dat');
fileID = fopen(filename, 'r');
pX = fscanf(fileID, '%f', [1, Inf]);

filename = strcat(folder, '/gridCellY.dat');
fileID = fopen(filename, 'r');
pY = fscanf(fileID, '%f', [1, Inf]);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

v = VideoWriter('video.avi','Motion JPEG AVI');
open(v);

fullname = strcat(folder, '/porosity');
fullname = strcat(fullname, '*.dat');
fstruct1 = dir(fullname);
fcell1 = struct2cell(fstruct1);

loops = numel(fstruct1)

h = figure;

set(gcf, 'Position',[50 50 1800 700]);

for k=1:loops

filename = strcat(folder, '/porosity');
filename = strcat(filename, string(k));
filename = strcat(filename, '.dat');
fileID = fopen(filename, 'r');
porosity = fscanf(fileID, '%f ', [1, Inf]);
porosity = reshape(porosity, M, N);

subplot(1,3,1)
plot(porosity(2,:), pY(2,:));
title('porosity')
ylabel('Depth');

filename = strcat(folder, '/darcypressure');
filename = strcat(filename, string(k));
filename = strcat(filename, '.dat');
fileID = fopen(filename, 'r');
darcypressure = fscanf(fileID, '%f ', [1, Inf]);
darcypressure = reshape(darcypressure, M, N);

subplot(1,3,2)
plot(darcypressure(2,:), pY(2,:));
title('darcypressure')
ylabel('Depth');

filename = strcat(folder, '/stokespressure');
filename = strcat(filename, string(k));
filename = strcat(filename, '.dat');
fileID = fopen(filename, 'r');
stokespressure = fscanf(fileID, '%f ', [1, Inf]);
stokespressure = reshape(stokespressure, M, N);

subplot(1,3,3)
plot(stokespressure(2,:), pY(2,:));
title('stokespressure')
ylabel('Depth');

pause
G = getframe(gcf);

writeVideo(v,G);

end

close(v);
