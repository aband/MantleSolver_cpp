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
fileID = fscanf(fileID);
porosity = fopen(fileID, '%f', [1,Inf]);


pause
G = getframe(gcf);

writeVideo(v,G);

end

close(v);
