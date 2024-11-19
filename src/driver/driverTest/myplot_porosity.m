%Plot porosity only
function [] = myplot_porosity(M, N)

fstruct = dir('build/*porosity*.dat');
fstruct = rmfield(fstruct,'folder');
fstruct = rmfield(fstruct,'date')
fstruct = rmfield(fstruct,'bytes')
fstruct = rmfield(fstruct,'isdir')
fstruct = rmfield(fstruct,'datenum')
fcell = struct2cell(fstruct);
for k=1:numel(fstruct)

fieldID = fopen(strcat('build/',fcell{k}), 'r')

end


fileID = fopen('build/gridCellX.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridCellY.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/porosity.dat','r');
poro = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/porosityFinal.dat','r');
finalporo = fscanf(fileID, '%f', [1, Inf]);

fclose(fileID);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);
poro = reshape(poro, M, N);
finalporo = reshape(finalporo, M, N);

co = [1 0 0; 0 1 0; 0 0 1];

figure
contourf(pX, pY, poro,20);
title("t = 0")
ylabel("depth (Dimensionless)");
set(gcf, 'Position',[0 0 300 600])
colorbar

figure
contourf(pX, pY, finalporo,20);
title("t = 0.2")
ylabel("depth (Dimensionless)")
set(gcf, 'Position',[0 0 300 600])
colorbar
