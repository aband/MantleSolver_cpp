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

for k=1:numel(fstruct1)

fileID = fopen(strcat('build/',fcell1{k}), 'r');
stokesx = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen(strcat('build/',fcell2{k}), 'r');
stokesy = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen(strcat('build/',fcell3{k}), 'r');
darcyx = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen(strcat('build/',fcell4{k}), 'r');
darcyy = fscanf(fileID, '%f', [1,Inf]);

stokesx = reshape(stokesx, M, N);
stokesy = reshape(stokesy, M, N);
darcyx = reshape(darcyx, M, N);
darcyy = reshape(darcyy, M, N);

figure
subplot(1,2,1)
quiver(pX, pY, stokesx, stokesy);
title("Stokes Velocity")

subplot(1,2,2)
quiver(pX, pY, darcyx, darcyy);
title("Darcy Velocity")

end

fclose(fileID);
