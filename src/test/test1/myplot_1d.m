function [] = myplot_1d(M, N)

% Read grid files 

fileID = fopen('build/gridCellX.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridCellY.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

v = VideoWriter('video.avi','Motion JPEG AVI');
open(v);


fstruct1 = dir('build/*porosity*.dat');
fcell1 = struct2cell(fstruct1);

fstruct2 = dir('build/*temperature*.dat');
fcell2 = struct2cell(fstruct2);

fstruct3 = dir('build/*stokesVy*.dat');
fcell3 = struct2cell(fstruct3);

fstruct4 = dir('build/*darcyVy*.dat');
fcell4 = struct2cell(fstruct4);

fstruct5 = dir('build/*stokesq*.dat');
fcell5 = struct2cell(fstruct5);

fstruct6 = dir('build/*darcyq*.dat');
fcell6 = struct2cell(fstruct6);

fstruct7 = dir('build/*CD*.dat');
fcell7 = struct2cell(fstruct7);

fstruct8 = dir('build/*HD*.dat');
fcell8 = struct2cell(fstruct8);

loops = numel(fstruct1)

h = figure;

for k=1:loops

% Porosity
filename = strcat('build/porosity',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
porosity = fscanf(fileID, '%f', [1,Inf]);
porosity = reshape(porosity, M, N);

subplot(1,8,1)
plot(porosity(2,:), pY(2,:))
title("Porosity")
ylabel("Depth");

% Temperature
filename = strcat('build/temperature',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
temp = fscanf(fileID, '%f', [1,Inf]);
temp = reshape(temp, M, N);

subplot(1,8,2)
plot(temp(2,:),pY(2,:));
title("Temperature");
%ylabel("Depth");

% Stokes Velocity in y direction
filename = strcat('build/stokesVy',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
stokesy = fscanf(fileID, '%f', [1,Inf]);
stokesy = reshape(stokesy, M, N);

subplot(1,8,3)
plot(stokesy(2,:),pY(2,:));
title("Stokes Velocity");
%ylabel("Depth");

% Unscaled darcy velocity in y driection
filename = strcat('build/darcyVy',string(k));
filename = strcat(filename,'.dat');
fileid = fopen(filename, 'r');
darcyy = fscanf(fileid, '%f', [1,inf]);
darcyy = reshape(darcyy, M, N);
unscaleddarcyy = darcyy.*porosity;

subplot(1,8,4)
plot(unscaleddarcyy(2,:),pY(2,:));
title("Darcy Velocity");

% Stokes pressure potential
filename = strcat('build/stokesq',string(k));
filename = strcat(filename,'.dat');
fileid = fopen(filename, 'r');
stokesq = fscanf(fileid, '%f', [1,inf]);
stokesq = reshape(stokesq, M, N);

subplot(1,8,5)
plot(stokesq(2,:),pY(2,:));
title("Stokes q");

% Darcy pressure potential
filename = strcat('build/darcyq',string(k));
filename = strcat(filename,'.dat');
fileid = fopen(filename, 'r');
darcyq = fscanf(fileid, '%f', [1,inf]);
darcyq = reshape(darcyq, M, N);

subplot(1,8,6)
plot(darcyq(2,:),pY(2,:));
title("Darcy q");

% dimensionless composition
filename = strcat('build/CD',string(k));
filename = strcat(filename,'.dat');
fileid = fopen(filename, 'r');
cd = fscanf(fileid, '%f', [1,inf]);
cd = reshape(cd, M, N);

subplot(1,8,7)
plot(cd(2,:),pY(2,:));
axis([-0.5 0.5 -0.2 0])
title("CD");

% dimensionless enthalpy
filename = strcat('build/HD',string(k));
filename = strcat(filename,'.dat');
fileid = fopen(filename, 'r');
hd = fscanf(fileid, '%f', [1,inf]);
hd = reshape(hd, M, N);

subplot(1,8,8)
plot(hd(2,:),pY(2,:));
title("HD");

set(gcf, 'Position',[50 50 1800 700]);

%pause
G = getframe(gcf);

writeVideo(v,G);
end

close(v);
