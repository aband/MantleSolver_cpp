function [] = myplot_pressure(M, N, L, start)

fileID = fopen('build/gridCellX.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridCellY.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

mid = floor(M/2) + 1;

v = VideoWriter('pressure.avi','Motion JPEG AVI');
open(v);

fstruct1 = dir('build/*porosity*.dat');
fcell1 = struct2cell(fstruct1);

loops = numel(fstruct1)

h = figure;

set(gcf, 'Position',[100 100 2010 1200])

mid = floor(M/2) + 1;

for kk=1:numel(fstruct1)

k = kk + start;

filename = strcat('build/qf', string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
data = fscanf(fileID, '%f', [1,Inf]);
data = reshape(data, M, N);

subplot(2,4,1)
surf(pX, pY, data);
title(filename)

subplot(2,4,1+4)
plot(data(mid,:),pY(mid,:));
ylim([-1*L,0.0])
title(filename);
ylabel("Depth");
xlabel("qf");

filename = strcat('build/qs', string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
data = fscanf(fileID, '%f', [1,Inf]);
data = reshape(data, M, N);

subplot(2,4,2)
surf(pX, pY, data);
title(filename)

subplot(2,4,2+4)
plot(data(mid,:),pY(mid,:));
ylim([-1*L,0.0])
title(filename);
ylabel("Depth");
xlabel("qs");

filename = strcat('build/rawstokesq', string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
data = fscanf(fileID, '%f', [1,Inf]);
data = reshape(data, M, N);

subplot(2,4,3)
surf(pX, pY, data);
title(filename)

subplot(2,4,3+4)
plot(data(mid,:),pY(mid,:));
ylim([-1*L,0.0])
title(filename);
ylabel("Depth");
xlabel("rawstokesq");

filename = strcat('build/rawdarcyq', string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
data = fscanf(fileID, '%f', [1,Inf]);
data = reshape(data, M, N);

subplot(2,4,4)
surf(pX, pY, data);
title(filename)

subplot(2,4,4+4)
plot(data(mid,:),pY(mid,:));
ylim([-1*L,0.0])
title(filename);
ylabel("Depth");
xlabel("rawdarcyq");

pause
F= getframe(gcf);

writeVideo(v,F);

end

fclose(fileID);
