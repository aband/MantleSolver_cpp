function [] = myplot_velocity(M,N)

% Read grid files
fileID = fopen('build/gridCellX.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridCellY.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

fstruct1 = dir('build/*stokesVx*.dat');
fstruct2 = dir('build/*stokesVy*.dat');
fstruct3 = dir('build/*darcyVx*.dat');
fstruct4 = dir('build/*darcyVy*.dat');

fstruct1 = rmfield(fstruct1,'folder');
fstruct1 = rmfield(fstruct1,'date');
fstruct1 = rmfield(fstruct1,'bytes');
fstruct1 = rmfield(fstruct1,'isdir');
fstruct1 = rmfield(fstruct1,'datenum');

fstruct2 = rmfield(fstruct2,'folder');
fstruct2 = rmfield(fstruct2,'date');
fstruct2 = rmfield(fstruct2,'bytes');
fstruct2 = rmfield(fstruct2,'isdir');
fstruct2 = rmfield(fstruct2,'datenum');

fstruct3 = rmfield(fstruct3,'folder');
fstruct3 = rmfield(fstruct3,'date');
fstruct3 = rmfield(fstruct3,'bytes');
fstruct3 = rmfield(fstruct3,'isdir');
fstruct3 = rmfield(fstruct3,'datenum');

fstruct4 = rmfield(fstruct4,'folder');
fstruct4 = rmfield(fstruct4,'date');
fstruct4 = rmfield(fstruct4,'bytes');
fstruct4 = rmfield(fstruct4,'isdir');
fstruct4 = rmfield(fstruct4,'datenum');

fcell1 = struct2cell(fstruct1);
fcell2 = struct2cell(fstruct2);
fcell3 = struct2cell(fstruct3);
fcell4 = struct2cell(fstruct4);

fstruct = dir('build/*porosity*.dat');
fstruct = rmfield(fstruct,'folder');
fstruct = rmfield(fstruct,'date');
fstruct = rmfield(fstruct,'bytes');
fstruct = rmfield(fstruct,'isdir');
fstruct = rmfield(fstruct,'datenum');
fcell = struct2cell(fstruct);

figure

for k=1:numel(fstruct1)

filename = strcat('build/stokesVx',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
stokesx = fscanf(fileID, '%f', [1,Inf]);

filename = strcat('build/stokesVy',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
stokesy = fscanf(fileID, '%f', [1,Inf]);

filename = strcat('build/darcyVx',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
darcyx = fscanf(fileID, '%f', [1,Inf]);

filename = strcat('build/darcyVy',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
darcyy = fscanf(fileID, '%f', [1,Inf]);

filename = strcat('build/porosity',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
poro   = fscanf(fileID, '%f', [1,Inf]);

stokesx = reshape(stokesx, M, N);
stokesy = reshape(stokesy, M, N);
darcyx = reshape(darcyx, M, N);
darcyy = reshape(darcyy, M, N);
poro   = reshape(poro,M,N);

subplot(1,3,1)
quiver(pX, pY, stokesx, stokesy);
title(["Stokes Velocity",num2str(k)])

subplot(1,3,2)
quiver(pX, pY, darcyx, darcyy);
title(["Scaled Darcy Velocity",num2str(k)])

subplot(1,3,3)
quiver(pX, pY, darcyx.*poro, darcyy.*poro);
title(["Unscaled Darcy Velocity",num2str(k)])

pause
F = getframe(gcf);
end

fclose(fileID);
