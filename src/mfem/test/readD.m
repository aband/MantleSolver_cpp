clc; clear

fileID = fopen('build/MatrixCheckAd.dat','r');

Ad = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/MatrixCheckBd.dat','r');

Bd = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/MatrixCheckgd1.dat','r');

gd1 = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/MatrixCheckgd2.dat','r');

gd2 = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

dofd1 = size(gd1,2);
dofd2 = size(gd2,2);

Ad = reshape(Ad,dofd1,dofd1);
Bd = reshape(Bd,dofd2,dofd1);

