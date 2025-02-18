% Plot phase transitions
function [] = myplot_phase(M, N)

% Input grid files
fileID = fopen('build/gridCellX.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridCellY.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

v = VideoWriter('phase.avi','Motion JPEG AVI');
open(v);

fstruct1 = dir('build/*porosity*.dat');
fcell1 = struct2cell(fstruct1);

fstruct2 = dir('build/*temperature*.dat');
fcell2 = struct2cell(fstruct2);

loops = numel(fstruct1)

h = figure;

for k=1:loops

filename = strcat('build/porosity',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
porosity = fscanf(fileID, '%f', [1,Inf]);

porosity = reshape(porosity, M, N);

subplot(1,8,1)
contourf(pX, pY, porosity,20);
title(strcat('porosity',string(k)));
ylabel("depth (Dimensionless)");
colorbar
%caxis([0,0.5])

% Add 1D plot
subplot(1,8,2)
plot(porosity(2,:),pY(2,:));
title("Porosity Distribution");
ylabel("Depth");
xlabel("Porosity");

filename = strcat('build/temperature',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
data = fscanf(fileID, '%f', [1,Inf]);

data = reshape(data, M, N);

subplot(1,8,3)
contourf(pX, pY, data,20);
title(strcat('temperature',string(k)))
ylabel("depth (Dimensionless)");
colorbar
%caxis([0,0.5])

% Add 1D plot
subplot(1,8,4)
plot(data(2,:),pY(2,:));
title("Temperature Distribution");
ylabel("Depth");
xlabel("Temperature");

%plot of HD
filename = strcat('build/HD',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
data = fscanf(fileID, '%f', [1,Inf]);

data = reshape(data, M, N);

subplot(1,8,5)
contourf(pX, pY, data,20);
title(strcat('HD',string(k)))
ylabel("depth (Dimensionless)");
colorbar
%caxis([0,0.5])

% Add 1D plot
subplot(1,8,6)
plot(data(2,:),pY(2,:));
title("HD Distribution");
ylabel("Depth");
xlabel("HD");

% Plot of CD
filename = strcat('build/CD',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
data = fscanf(fileID, '%f', [1,Inf]);

data = reshape(data, M, N);

subplot(1,8,7)
contourf(pX, pY, data,20);
title(strcat('CD',string(k)))
ylabel("depth (Dimensionless)");
colorbar
%caxis([0,0.5])

% Add 1D plot
subplot(1,8,8)
plot(data(2,:),pY(2,:));
title("CD Distribution");
ylabel("Depth");
xlabel("CD");

set(gcf, 'Position',[50 50 1800 700])

pause
F= getframe(gcf);

writeVideo(v,F);

end

fclose(fileID);
