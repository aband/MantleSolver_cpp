function [] = myplot_velocity(M,N)
% This plot generates beautiful pics for mid ocean ridges

fileID = fopen('build/gridCellX.dat','r');
cX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridCellY.dat','r');
cY = fscanf(fileID, '%f', [1,Inf]);

cX = reshape(cX, M, N);
cY = reshape(cY, M, N);

filename = strcat('build/porosity1.dat');
fileID = fopen(filename, 'r');
poro = fscanf(fileID, '%f', [1,Inf]);
poro = reshape(poro, M, N);


%figure
%set(gcf,'Position',[50 50 1000 500])
%x2 = linspace(min(min(cX)),max(max(cX)),200);
%y2 = linspace(min(min(cY)),max(max(cY)),200);
%[X2,Y2] = meshgrid(y2, x2);
%newporo = interp2(cX',cY',poro',X2,Y2);
%imagesc(X2,Y2,newporo')
%surf(cX,cY,poro,'LineStyle','none', 'FaceColor', 'interp');
%surf(cX,cY,poro,'LineStyle','none');
%colormap(turbo)
%view([0,90])
%colorbar
%title("Porosity Distrubition")
%xlabel('x')
%ylabel('y')

filename = strcat('build/qs1.dat');
fileID = fopen(filename, 'r');
qs = fscanf(fileID, '%f', [1,Inf]);
qs = reshape(qs, M, N);

filename = strcat('build/qf1.dat');
fileID = fopen(filename, 'r');
qf = fscanf(fileID, '%f', [1,Inf]);
qf = reshape(qf, M, N);

vM = 3*M;
vN = N+1;

% Read grid files
fileID = fopen('build/gaussgridx.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gaussgridy.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

% Extract grid
pX = reshape(pX, vM, vN);
pY = reshape(pY, vM, vN);

filename = strcat('build/darcyvely1.dat');
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, vM, vN);

filename = strcat('build/darcyvelx1.dat');
fileID = fopen(filename, 'r');
vx = fscanf(fileID, '%f', [1,Inf]);
vx = reshape(vx, vM, vN);

filename = strcat('build/stokesvely1.dat');
fileID = fopen(filename, 'r');
uy = fscanf(fileID, '%f', [1,Inf]);
uy = reshape(uy, vM, vN);

filename = strcat('build/stokesvelx1.dat');
fileID = fopen(filename, 'r');
ux = fscanf(fileID, '%f', [1,Inf]);
ux = reshape(ux, vM, vN);

filename = strcat('build/phasevy1.dat');
fileID = fopen(filename, 'r');
phasey = fscanf(fileID, '%f', [1,Inf]);
phasey = reshape(phasey, vM, vN);

filename = strcat('build/phasevx1.dat');
fileID = fopen(filename, 'r');
phasex = fscanf(fileID, '%f', [1,Inf]);
phasex = reshape(phasex, vM, vN);

pX  = pX(2:3:end, :);
pY  = pY(2:3:end, :);
vx  = vx(2:3:end, :);
vy  = vy(2:3:end, :);
ux  = ux(2:3:end, :);
uy  = uy(2:3:end, :);
phasex  = phasex(2:3:end, :);
phasey  = phasey(2:3:end, :);

figure
set(gcf,'Position',[50 50 1000 800])
%surf(cX,cY,qs,'LineStyle','none');
%hold on
quiver(pX, pY, ux, uy);

%[startx, starty] = meshgrid( 0.005,-0.08:0.004:-0.01)
%verts = stream2(pX',pY',ux',uy',startx,starty);
%streamline(verts)
%
%[startx, starty] = meshgrid(-0.005,-0.08:0.004:-0.01)
%verts = stream2(pX',pY',ux',uy',startx,starty);
%streamline(verts)
l = streamslice(pX',pY',ux',uy',1);
%l = streamslice(pX',pY',ux',uy');

set(l,'LineWidth',2);
set(l,'Color','k')

hold off
%colormap(turbo)
%view([0,90])
%colorbar
title("Solid Velocity")
xlabel('x')
ylabel('y')
xlim([-0.5,0.5])
%xlim([0.0,0.5])
ylim([-0.5,0.0])

%{
figure
set(gcf,'Position',[50 50 1000 500])
surf(cX,cY,qs,'LineStyle','none');
hold on
quiver(pX, pY, phasex, phasey);

%[startx, starty] = meshgrid( 0.005,-0.08:0.004:-0.01)
%verts = stream2(pX',pY',ux',uy',startx,starty);
%streamline(verts)
%
%[startx, starty] = meshgrid(-0.005,-0.08:0.004:-0.01)
%verts = stream2(pX',pY',ux',uy',startx,starty);
%streamline(verts)
l = streamslice(pX',pY',phasex',phasey',1.2);

set(l,'LineWidth',2);
set(l,'Color','k')

hold off
colormap(turbo)
view([0,90])
colorbar
title("Phase averaged velocity")
xlabel('x')
ylabel('y')
%xlim([-0.5,0.5])
xlim([0.0,0.5])
ylim([-0.5,0.0])

figure
set(gcf,'Position',[50 50 1000 500])
surf(cX,cY,qf,'LineStyle','none');
hold on
quiver(pX, pY, vx, vy);
l = streamslice(pX',pY',vx',vy',0.3);
set(l,'LineWidth',2);
set(l,'Color','k')
colormap(turbo)
view([0,90])
colorbar
title("Liquid Pressure Potential and Velocity")
xlabel('x')
ylabel('y')
%}
