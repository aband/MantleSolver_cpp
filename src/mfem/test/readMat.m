clc; clear

% Read matrix

fileID = fopen('build/MatrixCheckAs.dat','r');

As = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/MatrixCheckATest.dat','r');

ATest = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

fileID = fopen('build/MatrixCheckAd.dat','r');

Ad = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

fileID = fopen('build/MatrixCheckBs.dat','r');

Bs = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

fileID = fopen('build/MatrixCheckBd.dat','r');

Bd = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

fileID = fopen('build/MatrixCheckgs1.dat','r');

gs1 = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

fileID = fopen('build/MatrixCheckgs2.dat','r');

gs2 = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

fileID = fopen('build/MatrixCheckgd1.dat','r');

gd1 = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

fileID = fopen('build/MatrixCheckgd2.dat','r');

gd2 = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

fileID = fopen('build/MatrixCheckCs.dat','r');

Cs = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

fileID = fopen('build/MatrixCheckCd.dat','r');

Cd = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

fileID = fopen('build/MatCheckK.dat','r');

K = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);
% =====================================================================

dofs1 = size(gs1,2);
dofs2 = size(gs2,2);

dofd1 = size(gd1,2);
dofd2 = size(gd2,2);

As = reshape(As,dofs1,dofs1);
Bs = reshape(Bs,dofs2,dofs1);

ATest = reshape(ATest, 24, 24);

Ad = reshape(Ad,dofd1,dofd1);
Bd = reshape(Bd,dofd2,dofd1);

K = reshape(K,dofs2,dofs2);
Cs = reshape(Cs,dofs2,dofs2);
Cd = reshape(Cd,dofs2,dofs2);

A = [As, zeros(dofs1,dofd1);...
     zeros(dofd1,dofs1),Ad];

B = [Bs, zeros(dofs2,dofd1);...
     zeros(dofd2,dofs1), Bd];

C = [Cs, K; K, Cd];

F = [gs1';gd1'];

G = [gs2';gd2'];

tau1 = 10;
tau2 = 0.15;

tau = [tau1*eye(size(Cs)),zeros(dofs2,dofs2);
       zeros(dofs2,dofs2),tau2*eye(size(Cd))];
