clc; clear

% Read mesh points
% Read form gridX and gridY text file

fileID = fopen('build/gridX.txt','r');
X = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridY.txt','r');
Y = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);
