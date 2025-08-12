function [] = reconplot(M, N, mark, mytitle)

% Read grid files 

%fileID = fopen('savedrun/gridreconx.dat','r');
fileID = fopen('build/gridreconx.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

%fileID = fopen('savedrun/gridrecony.dat','r');
fileID = fopen('build/gridrecony.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, 3*M, 3*N);
pY = reshape(pY, 3*M, 3*N);

%filename = strcat('savedrun/reconSol',string(mark));
filename = strcat('build/reconSol',string(mark));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
sol = fscanf(fileID, '%f', [1,Inf]);
sol = reshape(sol, 3*M, 3*N);

figure
s = surf(pX, pY, sol)
%s.EdgeColor = 'none';
title(mytitle)
ylabel("y")
xlabel("x")
colormap(turbo)
%colorbar
axis equal
%xlim([0,3])
%ylim([0,1])
%zlim([0,1])
set(gcf, 'Position', [50 50 1000 800]); % Set position and size of the current figure

figure
mesh(pX,pY)

fig = figure 
ax = axes;
colormap(turbo)
caxis([0,1])
c = colorbar(ax)
ax.Visible = 'off'

exportgraphics(fig, 'colorbar.png','resolution',300)
