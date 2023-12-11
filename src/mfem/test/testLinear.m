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

% Enlarge system M

M = [M,zeros(dof1+dof2,1);zeros(1,dof1+dof2),1];

g = [g;0];

% ============== Checking ==========================

rank(B)
rank([A,B'])
