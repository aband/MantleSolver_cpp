clc; clear

fileID = fopen('build/gridCD.dat','r');
CD = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/gridHD.dat','r');
HD = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/TD.dat','r');
TD = fscanf(fileID, '%f', [1,Inf]);

fclose(fileID);

CD = reshape(CD, 51, 51);
HD = reshape(HD, 51, 51);
TD = reshape(TD, 51, 51);

% Surface plot of dimensionless temperature
% regarding dimensionless composition and dimensionless enthalpy
figure
surf(CD, HD, TD);

% ============================================================

% Plotting slice of (CD, HD, TD)



