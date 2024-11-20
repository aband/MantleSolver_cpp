% Plot porosity only
function [] = myplot_porosity(M, N)

% Input grid files
fileID = fopen('build/gridCellX.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridCellY.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

fstruct = dir('build/*porosity*.dat');
fstruct = rmfield(fstruct,'folder');
fstruct = rmfield(fstruct,'date');
fstruct = rmfield(fstruct,'bytes');
fstruct = rmfield(fstruct,'isdir');
fstruct = rmfield(fstruct,'datenum');
fcell = struct2cell(fstruct);

for k=1:numel(fstruct)

fileID = fopen(strcat('build/',fcell{k}), 'r');
data = fscanf(fileID, '%f', [1,Inf]);

data = reshape(data, M, N);

figure
contourf(pX, pY, data,20);
title(fcell{k})
ylabel("depth (Dimensionless)");
set(gcf, 'Position',[0 0 300 600])
colorbar

end

fclose(fileID);
