function [] = velplot(M, N, mark, title)

% Read grid files
fileID  = fopen('test/build/vertgaussgridx.dat');
vertpx  = fscanf(fileID, '%f', [1,Inf]);
fclose(fileID);

fileID  = fopen('test/build/vertgaussgridy.dat');
vertpy  = fscanf(fileID, '%f', [1,Inf]);
fclose(fileID);

fileID  = fopen('test/build/horigaussgridx.dat');
horipx  = fscanf(fileID, '%f', [1,Inf]);
fclose(fileID);

fileID  = fopen('test/build/horigaussgridy.dat');
horipy  = fscanf(fileID, '%f', [1,Inf]);
fclose(fileID);

% ===============================================

MM = M*3;
NN = N*3;

vertpx = reshape(vertpx, );
vertpy = reshape(vertpy, );
