function [] = myplot(M, N)

% Read grid files
fileID = fopen('build/gridCellX.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridCellY.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/porosity.dat','r');
poro = fscanf(fileID, '%f', [1,Inf]);

% Read velocity files(unscaled yet)
fileID = fopen('build/stokesVx.dat','r');
stokesvx = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/stokesVy.dat','r');
stokesvy = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/darcyVx.dat','r');
darcyvx = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/darcyVy.dat','r');
darcyvy = fscanf(fileID, '%f', [1,Inf]);

%Read pressure files(unscaled yet)
fileID = fopen('build/stokesp.dat','r');
stokesp = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/darcyp.dat','r');
darcyp = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);
poro = reshape(poro, M, N);

stokesvx = reshape(stokesvx, M, N);
stokesvy = reshape(stokesvy, M, N);
darcyvx = reshape(darcyvx, M, N);
darcyvy = reshape(darcyvy, M, N);

darcyp = reshape(darcyp, M, N);
stokesp = reshape(stokesp, M, N);

figure
contour(pX, pY, poro,20);
title("Porosity Contour")

figure 
quiver(pX, pY, stokesvx, stokesvy);
title('Stokes Velocity');

figure
quiver(pX, pY, darcyvx, darcyvy);
title('Darcy Velocity');

figure 
surf(pX, pY, darcyp);
title('Unscaled darcy pressure');

figure
surf(pX, pY, stokesp);
title('Stokes pressure');
