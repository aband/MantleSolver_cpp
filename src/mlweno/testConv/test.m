clc; clear

X = readmatrix("build/gridX.txt");
Y = readmatrix("build/gridY.txt");

X = (X(:,1:end-1) + X(:,2:end))/2;
Y = (Y(1:end-1,:) + Y(2:end,:))/2;

X = X(1:end-1,:);
Y = Y(:,1:end-1);

t = 0;

M = readmatrix("build/initial.txt");

figure
h = surf(X,Y,M);
set(h,'linestyle','none')
title(['t = ', num2str(t), ' s']);
