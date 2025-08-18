function [] = reconplot(M, N, mark, mytitle)

% Read grid files

fileID = fopen('build/exactgridx.dat');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/exactgridy.dat');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, 3*M, 3*N);
pY = reshape(pY, 3*M, 3*N);

%{
filename = strcat('build/exactsol',string(mark));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
sol = fscanf(fileID, '%f', [1,Inf]);
sol = reshape(sol, 3*M, 3*N);

figure
s = surf(pX, pY, sol)
s.EdgeColor = 'none';
title('exact')
ylabel("y")
xlabel("x")
colormap(turbo)
colorbar
caxis([0,1])
%}

filename = strcat('build/reconSol',string(mark));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
sol2 = fscanf(fileID, '%f', [1,Inf]);
sol2 = reshape(sol2, 3*M, 3*N);

figure
s = surf(pX, pY, sol2)
s.EdgeColor = 'none';
view(360,0)
title(mytitle)
ylabel("y")
xlabel("x")
colormap(turbo)
set(gcf, 'Position',[50 50 350 1200]);
ylim([0,1.2])
xlim([0,0.23])
zlim([-0.2,1.2])
axis([0 0.23 0 1.2])
%axis equal

%colorbar
%caxis([0,1])

fig = figure 
ax = axes;
colormap(turbo)
caxis([0,1])
c = colorbar(ax)
ax.Visible = 'off'

exportgraphics(fig, 'colorbar.png','resolution',300)
