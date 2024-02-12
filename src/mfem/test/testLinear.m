% Test generated linear system with matlab first

clc; clear
% First read everything from plain text
fileID = fopen('build/MatrixCheckA.dat','r');

A = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

fileID = fopen('build/MatrixCheckB.dat','r');

B = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

fileID = fopen('build/MatrixCheckg1.dat','r');

g1 = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

fileID = fopen('build/MatrixCheckg2.dat','r');

g2 = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

fileID = fopen('build/VecCheckg.dat','r');

og = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

fileID = fopen('build/MatrixCheckKg.dat','r');

Kg = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

%fileID = fopen('build/MatrixCheckFullA.dat','r');
%fullA = fscanf(fileID, '%f', [1,Inf]);
%fclose(fileID);

% ==================================================

dof1 = size(g1,2);
dof2 = size(g2,2);
%dof3 = size(og,2);

A = reshape(A,dof1,dof1);
B = reshape(B,dof2,dof1);
%Kg = reshape(Kg,dof1,dof1);

g = [g1';g2'];

Z = zeros(dof2,dof2);

M = [A,B';B,Z];

% ==================================================
% Test 1
% Enlarge system M

%lg = [zeros(dof1,1);ones(dof2,1)];

%Ml = [M,lg;lg',0];

%gl = [g;0];

%Ml\gl;

% ==================================================
% Test 2
% Uzawa iteration
r = 1.0;
MaxIter = 3000;
iter = 0;

F = g1';
G = g2';

x = zeros(dof1,1);
y = zeros(dof2,1);

while (r > 1e-8) && (iter < MaxIter)

    iter + 1;

    tmp1 = A\(F - (A*x - B'*y));

    x = x + tmp1;

    tmp2 = -B*x+G;
  
    y = y + 1*tmp2;

    iter = iter +1;

    r = norm(tmp1) + norm(tmp2);
end

% ======== Test of =========

%ux = @(x,y) -x./(x.^2+y.^2);
%uy = @(x,y) -y./(x.^2+y.^2);
ux = @(x,y)  x;
uy = @(x,y) -y;

cx = linspace(-1,1,20);
cy = linspace(-1,1,20);

[X,Y] = meshgrid(cx,cy);
U = ux(X,Y);
V = uy(X,Y);

%quiver(X,Y,U,V);
