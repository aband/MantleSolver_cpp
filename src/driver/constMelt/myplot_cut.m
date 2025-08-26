function [] = myplot_cut(M,N, cut)

% Print porosity and velocity at the given time stamp

fileID = fopen('build/gridCellX.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridCellY.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

filename = strcat('build/porosity', string(cut));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
poro = fscanf(fileID, '%f', [1,Inf]);
poro = reshape(poro, M, N);

filename = strcat('build/stokesvx',string(cut));
filename = strcat(filename,'.dat');
fileid = fopen(filename, 'r');
stokesx = fscanf(fileid, '%f', [1,inf]);

filename = strcat('build/stokesVy',string(cut));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
stokesy = fscanf(fileID, '%f', [1,Inf]);

filename = strcat('build/darcyVx',string(cut));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
darcyx = fscanf(fileID, '%f', [1,Inf]);

filename = strcat('build/darcyVy',string(cut));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
darcyy = fscanf(fileID, '%f', [1,Inf]);

stokesy = reshape(stokesy, M, N);
stokesx = reshape(stokesx, M, N);
darcyy = reshape(darcyy, M, N);
darcyx = reshape(darcyx, M, N);

unscaleddarcyy = darcyy.*poro;

% ============================================================
subplot(1,3,2)
plot(poro(2,:), pY(2,:));
title(filename)
ylabel("Depth");
axis([])

subplot(1,3,2)
quiver(pX, pY, darcyx.*poro, unscaleddarcyy);
title(["Darcy Velocity"])


%axis equal
