function [] = myplot_velocity(M,N)
% This function generates beautiful pics

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

filename = 'build/exactql';
fileID = fopen(filename, 'r');
exactql = fscanf(fileID, '%f', [1,Inf]);
exactql = reshape(exactql, M, N);

filename = 'build/exactqs';
fileID = fopen(filename, 'r');
exactqs = fscanf(fileID, '%f', [1,Inf]);
exactqs = reshape(exactqs, M, N);

filename = strcat('build/qs1.dat');
fileID = fopen(filename, 'r');
qs = fscanf(fileID, '%f', [1,Inf]);
qs = reshape(qs, M, N);

filename = strcat('build/qf1.dat');
fileID = fopen(filename, 'r');
qf = fscanf(fileID, '%f', [1,Inf]);
qf = reshape(qf, M, N);

%exactql(2,1:N/2) = exactql(2,1:N/2) + mean(qs(1,:)); 
%exactql(2,N/2:N) = exactql(2,N/2:N) + mean(qs(1,:)); 
exactql(2,:) = exactql(2,:) - mean(qf(1,:)); 
exactqs(2,:) = exactqs(2,:) + mean(qs(1,:));
area = 4/N * 0.2

figure
set(gcf, 'Position',[50 50 250 700]);
hold on
plot(exactqs(2,:), cY(2,:), '-c', 'LineWidth', 1);
plot(exactql(2,:), cY(2,:), '-k', 'LineWidth', 1);
plot(qs(1,:), cY(1,:),'--r','LineWidth',2);
plot(-qf(1,:), cY(1,:),'--b','LineWidth',2);
hold off

sqrt(sum((abs(exactqs(2,:) - qs(2,:))).^2))*area
sqrt(sum((abs(exactql(2,:) + qf(2,:))).^2))*area

legend("exact q_s", "exact q_l", "q_s", "q_l");

figure
set(gcf, 'Position',[50 50 250 700]);
hold on
plot(poro(2,1:N/2+1), cY(2,1:N/2+1), '-k', 'LineWidth', 2);

plot(poro(2,N/2+1:N), cY(2,N/2+1:N), '-k', 'LineWidth', 2);
hold off
xlim([-0.002, 0.042])
%xlim([-0.0002, 0.0042])

% Read grid files
fileID = fopen('build/gaussgridx.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gaussgridy.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

M = 3*M;
N = N+1;

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

filename = ('build/exactv.dat');
fileID = fopen(filename,'r')
exact = fscanf(fileID, '%f', [1,Inf]);
exact = reshape(exact, M, N);

filename = strcat('build/darcyvely1.dat');
fileID = fopen(filename, 'r');
vy = fscanf(fileID, '%f', [1,Inf]);
vy = reshape(vy, M, N);

filename = strcat('build/darcyvelx1.dat');
fileID = fopen(filename, 'r');
vx = fscanf(fileID, '%f', [1,Inf]);
vx = reshape(vx, M, N);

filename = strcat('build/stokesvely1.dat');
fileID = fopen(filename, 'r');
uy = fscanf(fileID, '%f', [1,Inf]);
uy = reshape(uy, M, N);

filename = strcat('build/stokesvelx1.dat');
fileID = fopen(filename, 'r');
ux = fscanf(fileID, '%f', [1,Inf]);
ux = reshape(ux, M, N);

% ==========================================

figure
set(gcf, 'Position',[50 50 250 700]);

hold on
plot(exact(4,:), pY(4,:), '-k', 'LineWidth',1 );
plot(-1*exact(4,:), pY(4,:),'-k' ,'LineWidth',1);
plot(-1*uy(4,:), pY(4,:),'--r','LineWidth',2);
plot(vy(4,:), pY(4,:),'--b','LIneWidth',2);
hold off

legend("exact Stokes vel", "exact Darcy vel", "Stokes vel v", "Darcy vel u",'Location','north');

% ==========================================

figure
set(gcf, 'Position',[50 50 250 700]);

subplot(2,1,1)
quiver(pX, pY, vx, vy);
title(["darcy velocity u"])

subplot(2,1,2)
quiver(pX, pY, ux, uy);
title(["Stokes velocity v"])

end
