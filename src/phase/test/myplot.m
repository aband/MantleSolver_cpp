clc; clear

fileID = fopen('build/gridCD.dat','r');
CD = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/gridHD.dat','r');
HD = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/TD.dat','r');
TD = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/Vf.dat','r');
Vf = fscanf(fileID, '%f', [1,Inf]);

seed = 101;
fclose(fileID);

CD = reshape(CD, seed, seed);
HD = reshape(HD, seed, seed);
TD = reshape(TD, seed, seed);
Vf = reshape(Vf, seed, seed);

figure
h = surf(CD, HD, TD);
get(h)
%set(h,'linestyle','none','facecolor',[0 0.4470 0.7410]);
set(h,'linestyle','none','facecolor','interp');
light("Style","local","Position",[0 0 10]);
title("Composition-Enthalpy-Temperature");
xlabel("Composition");
ylabel("Enthalpy");
zlabel("Temperature");

figure
g = surf(CD, HD, Vf);
set(g,'linestyle','none','facecolor','interp');
light("Style","local","Position",[0 0 10]);
title("Composition-Enthalpy-VolumeFraction");
xlabel("Composition");
ylabel("Enthalpy");
zlabel("VolumrFraction");

% ============================================================
