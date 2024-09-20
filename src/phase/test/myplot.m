clc; clear

fileID = fopen('build/gridCD.dat','r');
CD = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/gridHD.dat','r');
HD = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/TD.dat','r');
TD = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/Vf.dat','r');
Vf = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/TDz.dat','r');
TDz = fscanf(fileID, '%f', [1,Inf]);
fclose(fileID);

CD = reshape(CD, 51, 51);
HD = reshape(HD, 51, 51);
TD = reshape(TD, 51, 51);
Vf = reshape(Vf, 51, 51);

% Surface plot of dimensionless temperature
% regarding dimensionless composition and dimensionless enthalpy
figure
surf(CD, HD, TD);
title("Composition-Enthalpy-Temperature");
figure
surf(CD, HD, Vf);
title("Composition-Enthalpy-VolumeFraction");
% ============================================================

gamma = 10^(-7);
rho   = 3000;
g     = 10;
T1    = 2053;
Te    = 1480;

param = gamma * rho * g / (T1-Te);

z = linspace(0,60000,100);
figure
plot(z, z*param,z,TDz)
hold on
xlabel('z(m)')
ylabel('TD')
legend("Eutectic point","Actual Temp")
% Plotting slice of (CD, HD, TD)


