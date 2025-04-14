function [] = myplot_pressure(M, N, L, start)

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

loops = numel(fstruct1)

h = figure;

set(gcf, 'Position',[100 100 1210 693])

for kk=1:numel(fstruct1)

k = kk + start

fileID = fopen(strcat('build/',fcell1{k}), 'r');
stokesq = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen(strcat('build/',fcell2{k}), 'r');
darcyq = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen(strcat('build/',fcell3{k}), 'r');
stokesp = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen(strcat('build/',fcell4{k}), 'r');
darcyp = fscanf(fileID, '%f', [1,Inf]);

stokesq = reshape(stokesq, M, N);
darcyq = reshape(darcyq, M, N);
stokesp = reshape(stokesp, M, N);
darcyp = reshape(darcyp, M, N);

figure
subplot(2,2,1)
surf(pX, pY, stokesq);
title(["Stokes Pressure Potential",num2str(k)])

subplot(2,2,2)
surf(pX, pY, darcyq);
title(["Darcy Pressure Potential",num2str(k)])

subplot(2,2,3)
surf(pX, pY, stokesp);
title(["Stokes Pressure",num2str(k)])

subplot(2,2,4)
surf(pX, pY, darcyp);
title(["Darcy Pressure",num2str(k)])

figure
subplot(2,2,1)
plot(stokesq(1,:),pY(2,:));
title(["Stokes Pressure Potential",num2str(k)])
ylabel("Depth");

subplot(2,2,2)
plot(darcyq(1,:),pY(2,:));
title(["Darcy Pressure Potential",num2str(k)])
ylabel("Depth")

subplot(2,2,3)
plot(stokesp(1,:), pY(2,:));
title(["Stokes Pressure",num2str(k)])
ylabel("Depth")

subplot(2,2,4)
plot(darcyp(1,:),pY(2,:));
title(["Darcy Pressure",num2str(k)])
ylabel("Depth")

end

fclose(fileID);
