function [] = shape(N)

fileID = fopen('build/gridX.dat','r');
X = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridY.dat','r');
Y = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/R1.dat','r');
R1 = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/R.dat','r');
R = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/phie.dat','r');
phie = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/phiv.dat','r');
phiv = fscanf(fileID, '%f', [1,Inf]);


X = reshape(X,N,N);
Y = reshape(Y,N,N);
R1 = reshape(R1,N,N);
R = reshape(R,N,N);
phie = reshape(phie,N,N);
phiv = reshape(phiv,N,N);



%figure
%surf(X,Y,R1);

figure
surf(X,Y,R,'LineStyle','none');
hold on
plot([0.1 0.8 1.2 -0.05 0.1], [-0.2 0.1 0.95 1.03, -0.2], 'k-');
hold off
title("$$R_0$$",'interpreter','latex')
view([10,30]);

figure
surf(X,Y,R,'LineStyle','none');
hold on
plot([0.1 0.8 1.2 -0.05 0.1], [-0.2 0.1 0.95 1.03, -0.2], 'k-');
hold off
title("$$R_0$$",'interpreter','latex')
view([0,90]);

figure
surf(X,Y,phie,'LineStyle','none');
hold on
plot([0.1 0.8 1.2 -0.05 0.1], [-0.2 0.1 0.95 1.03, -0.2], 'k-');
hold off
title("$$\varphi_{e,0}$$",'interpreter','latex')
view([10,30]);

figure
surf(X,Y,phie,'LineStyle','none');
hold on
plot([0.1 0.8 1.2 -0.05 0.1], [-0.2 0.1 0.95 1.03, -0.2], 'k-');
hold off
title("$$\varphi_{e,0}$$",'interpreter','latex')
view([0,90]);

figure
surf(X,Y,phiv,'LineStyle','none');
hold on
plot([0.1 0.8 1.2 -0.05 0.1], [-0.2 0.1 0.95 1.03, -0.2], 'k-');
hold off
title("$$\varphi_{v,0}$$",'interpreter','latex')
view([10,30]);

figure
surf(X,Y,phiv,'LineStyle','none');
hold on
plot([0.1 0.8 1.2 -0.05 0.1], [-0.2 0.1 0.95 1.03, -0.2], 'k-');
hold off
title("$$\varphi_{v,0}$$",'interpreter','latex')
view([0,90]);


