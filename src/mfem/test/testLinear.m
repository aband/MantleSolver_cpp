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

% ==================================================

dof1 = size(g1,2);
dof2 = size(g2,2);

A = reshape(A,dof1,dof1);
B = reshape(B,dof2,dof1);

g = [g1';g2'];

M = [A,B';B,zeros(dof2,dof2)];

% ==================================================
% Test 1
% Enlarge system M

lg = [zeros(dof1,1);ones(dof2,1)];

Ml = [M,lg;lg',0];

gl = [g;0];

Ml\gl;

% ==================================================
% Test 2
% Uzawa iteration

z = A\g1';

S = B/A * B'; % Schur complement

y = S\(B*z);

x = A\(z - B'*y);

% ===================================================
% Test 3
% Preconditioned MINRES method


