function [] = shape(N)

fileID = fopen('build/gridX1.dat','r');
X1 = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridY1.dat','r');
Y1 = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/valx1.dat','r');
vX1 = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/valy1.dat','r');
vY1 = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/ualx1.dat','r');
uX1 = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/ualy1.dat','r');
uY1 = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/sc1.dat','r');
sc1 = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/scc1.dat','r');
scc1 = fscanf(fileID, '%f', [1,Inf]);

X1 = reshape(X1,N,N);
Y1 = reshape(Y1,N,N);
vX1 = reshape(vX1,N,N);
vY1 = reshape(vY1,N,N);
uX1 = reshape(uX1,N,N);
uY1 = reshape(uY1,N,N);

sc1 = reshape(sc1,N,N);
scc1 = reshape(scc1,N,N);

fileID = fopen('build/gridX2.dat','r');
X2 = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridY2.dat','r');
Y2 = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/valx2.dat','r');
vX2 = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/valy2.dat','r');
vY2 = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/ualx2.dat','r');
uX2 = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/ualy2.dat','r');
uY2 = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/sc2.dat','r');
sc2 = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/scc2.dat','r');
scc2 = fscanf(fileID, '%f', [1,Inf]);

X2 = reshape(X2,N,N);
Y2 = reshape(Y2,N,N);
vX2 = reshape(vX2,N,N);
vY2 = reshape(vY2,N,N);
uX2 = reshape(uX2,N,N);
uY2 = reshape(uY2,N,N);

sc2 = reshape(sc2,N,N);
scc2 = reshape(scc2,N,N);

figure
quiver(X1,Y1,vX1,vY1);
hold on
plot([X1(1,1) X1(N,1) X1(N,N) X1(1,N) X1(1,1)], [Y1(1,1) Y1(N,1) Y1(N,N) Y1(1,N) Y1(1,1)], 'k-');
plot([X2(1,1) X2(N,1) X2(N,N) X2(1,N) X2(1,1)], [Y2(1,1) Y2(N,1) Y2(N,N) Y2(1,N) Y2(1,1)], 'k-');
quiver(X2,Y2,vX2,vY2);
hold off
title("$$\varphi_{l}$$",'interpreter','latex');
xlabel("x");
ylabel("y");

figure
surf(X1,Y1,sc1,'LineStyle','none','FaceColor','interp');
hold on
surf(X2,Y2,sc2,'LineStyle','none','FaceColor','interp');
hold off
colormap turbo
xlabel("x");
ylabel("y");

figure
quiver(X1,Y1,uX1,uY1);
hold on
plot([X1(1,1) X1(N,1) X1(N,N) X1(1,N) X1(1,1)], [Y1(1,1) Y1(N,1) Y1(N,N) Y1(1,N) Y1(1,1)], 'k-');
plot([X2(1,1) X2(N,1) X2(N,N) X2(1,N) X2(1,1)], [Y2(1,1) Y2(N,1) Y2(N,N) Y2(1,N) Y2(1,1)], 'k-');
quiver(X2,Y2,uX2,uY2);
hold off
title("$$\varphi_{c}$$",'interpreter','latex');
xlabel("x");
ylabel("y");

figure
surf(X1,Y1,scc1,'LineStyle','none','FaceColor','interp');
hold on
surf(X2,Y2,scc2,'LineStyle','none','FaceColor','interp');
hold off
colormap turbo
xlabel("x");
ylabel("y");

