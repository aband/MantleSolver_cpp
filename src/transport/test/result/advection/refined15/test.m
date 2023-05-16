clc; clear

X = readmatrix("gridX.txt");
Y = readmatrix("gridY.txt");

X = (X(:,1:end-1) + X(:,2:end))/2;
Y = (Y(1:end-1,:) + Y(2:end,:))/2;

X = X(1:end-1,:);
Y = Y(:,1:end-1);

t = 0;

M = readmatrix("initial.txt");
N = readmatrix("final.txt");

figure
h = surf(X,Y,M);
set(h,'linestyle','none')
title(['t = ', num2str(t), ' s']);

figure
h = surf(X,Y,N);
set(h,'linestyle','none')
title(['t = ', num2str(1.5), ' s']);
