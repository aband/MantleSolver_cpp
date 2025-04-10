% Plot phase transitions
function [] = myplot_phase(M, N, L)

% Input grid files
fileID = fopen('build/gridCellX.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridCellY.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

mid = floor(M/2) + 1;

v = VideoWriter('phase.avi','Motion JPEG AVI');
open(v);

fstruct1 = dir('build/*porosity*.dat');
fcell1 = struct2cell(fstruct1);

fstruct2 = dir('build/*temperature*.dat');
fcell2 = struct2cell(fstruct2);

loops = numel(fstruct1)

h = figure;

set(gcf, 'Position',[100 100 1210 693])

for k=1:loops

filename = strcat('build/porosity',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
porosity = fscanf(fileID, '%f', [1,Inf]);
porosity = reshape(porosity, M, N);

%subplot(1,5,1)
%contourf(pX, pY, porosity,20);
%title(strcat('porosity',string(k)));
%ylabel("depth (Dimensionless)");
%colorbar
%caxis([0,0.5])

% Add 1D plot
subplot(1,6,1)
plot(porosity(mid,:),pY(mid,:));
title("Porosity");
ylabel("Depth");
xlabel("Porosity");

filename = strcat('build/temperature',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
data = fscanf(fileID, '%f', [1,Inf]);
data = reshape(data, M, N);

filename = strcat('build/meltT',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
mT = fscanf(fileID, '%f', [1,Inf]);
mT = reshape(mT, M, N);

%subplot(1,8,3)
%contourf(pX, pY, data,20);
%title(strcat('temperature',string(k)))
%ylabel("depth (Dimensionless)");
%colorbar
%caxis([0,0.5])

% Add 1D plot
subplot(1,6,2)
plot(data(mid,:),pY(mid,:));
hold on 
plot(mT(mid,:),pY(mid,:), 'o');
hold off
title("Temperature");
ylabel("Depth");
xlabel("Temperature");

%plot of HD
filename = strcat('build/HD',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
data = fscanf(fileID, '%f', [1,Inf]);
data = reshape(data, M, N);

%subplot(1,8,5)
%contourf(pX, pY, data,20);
%title(strcat('HD',string(k)))
%ylabel("depth (Dimensionless)");
%colorbar
%caxis([0,0.5])

% Add 1D plot
subplot(1,6,3)
plot(data(mid,:),pY(mid,:));
title("HD Distribution");
ylabel("Depth");
xlabel("HD");

% Plot of CD
filename = strcat('build/CD',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
data = fscanf(fileID, '%f', [1,Inf]);
data = reshape(data, M, N);

%subplot(1,8,7)
%contourf(pX, pY, data,20);
%title(strcat('CD',string(k)))
%ylabel("depth (Dimensionless)");
%colorbar
%caxis([0,0.5])

% Add 1D plot
subplot(1,6,4)
plot(data(mid,:),pY(mid,:));
axis([0.00,0.11, -1*L, 0.0])
title("CD Distribution");
ylabel("Depth");
xlabel("CD");

% Plot of phase split region
filename = strcat('build/phase', string(k));
filename = strcat(filename, '.dat');
fileID   = fopen(filename, 'r');
data     = fscanf(fileID, '%f', [1, Inf]);
data = reshape(data, M, N);

subplot(1,6,5);
%contourf(pX, pY, data, 3);
stairs(data(mid,:), pY(mid,:));
axis([0, 4, -1*L,0.0])
title(strcat('Phase split'));

% Plot of phase split region
filename = strcat('build/opx', string(k));
filename = strcat(filename, '.dat');
fileID   = fopen(filename, 'r');
data     = fscanf(fileID, '%f', [1, Inf]);
data = reshape(data, M, N);

subplot(1,6,6);
plot(data(mid,:),pY(mid,:));
ylim([-1*L,0.0])
xlim([0.00, 0.12])
title(filename);
ylabel("Depth");
xlabel("opx");

pause
F= getframe(gcf);

writeVideo(v,F);

end

fclose(fileID);
