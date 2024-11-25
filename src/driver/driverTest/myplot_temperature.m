function [] = myplot_temperature(M, N)

% Input grid files
fileID = fopen('build/gridCellX.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridCellY.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

fstruct = dir('build/*temperature*.dat');
fstruct = rmfield(fstruct,'folder');
fstruct = rmfield(fstruct,'date');
fstruct = rmfield(fstruct,'bytes');
fstruct = rmfield(fstruct,'isdir');
fstruct = rmfield(fstruct,'datenum');
fcell = struct2cell(fstruct);

figure

for k=1:numel(fstruct)

filename = strcat('build/temperature',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
data = fscanf(fileID, '%f', [1,Inf]);

data = reshape(data, M, N);

subplot(1,2,1)
contourf(pX, pY, data,20);
title(strcat('temperature',string(k)))
ylabel("depth (Dimensionless)");
set(gcf, 'Position',[50 50 400 700])
colorbar
caxis([0,0.5])

% Add 1D plot
subplot(1,2,2)
plot(data(2,:),pY(2,:));
title("Temperature Distribution");
ylabel("Depth");
xlabel("Temperature");

pause
F = getframe(gcf);

end

fclose(fileID);
