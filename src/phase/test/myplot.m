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
fileID = fopen('build/Vfz.dat','r');
Vfz = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/TDz2.dat','r');
TDz2 = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/Vfz2.dat','r');
Vfz2 = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/HDzgrid.dat','r');
HDzgrid = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/Pzgrid.dat','r');
Pzgrid = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/TDz3.dat','r');
TDz3 = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/Vfz3.dat','r');
Vfz3 = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/TDz4.dat','r');
TDz4 = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/Vfz4.dat','r');
Vfz4 = fscanf(fileID, '%f', [1,Inf]);


fclose(fileID);

CD = reshape(CD, 51, 51);
HD = reshape(HD, 51, 51);
TD = reshape(TD, 51, 51);
Vf = reshape(Vf, 51, 51);

TDz2 = reshape(TDz2, 100, 100);
Vfz2 = reshape(Vfz2, 100, 100);

HDzgrid = reshape(HDzgrid, 100, 100);
Pzgrid = reshape(Pzgrid, 100, 100);

% Surface plot of dimensionless temperature
% regarding dimensionless composition and dimensionless enthalpy
%figure
%surf(CD, HD, TD);
%title("Composition-Enthalpy-Temperature");
%figure
%surf(CD, HD, Vf);
%title("Composition-Enthalpy-VolumeFraction");

% ============================================================

gamma = 10^(-7);
rho   = 3000;
g     = 10;
T1    = 2053;
Te    = 1480;

param = gamma * rho * g / (T1-Te);

z = linspace(0,60000,100);
figure
plot(z, z*param,'linewidth',1.5)
hold on
plot(z, TDz,'linewidth',1.5);
plot(z, TDz3,'linewidth',1.5);
plot(z, TDz4,'linewidth',1.5);

xlabel('z(m)')
ylabel('TD')
legend("Eutectic point","CD = 0.2","CD = 0.4", "CD = 0.6")
camroll(270)
set(gcf,'Position',[0 0 300 600])
hold off

figure
plot(z,Vfz,'linewidth',1.5)
hold on
plot(z,Vfz3,'linewidth',1.5);
plot(z,Vfz4,'linewidth',1.5)
xlabel('z(m)')
ylabel('Vf')
legend("CD = 0.2", "CD = 0.4", "CD = 0.6",'Location','southeast');
camroll(270)
set(gcf,'Position',[0 0 300 600])
hold off
% ============================================================

%figure
%surf(HDzgrid, Pzgrid, TDz2);

%figure
%surf(HDzgrid, Pzgrid, Vfz2);


