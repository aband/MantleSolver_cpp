clc; clear

fileID = fopen('build/gridCD.dat','r');
CD = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/gridHD.dat','r');
HD = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/TD.dat','r');
TD = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/Vf.dat','r');
Vf = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/dTdH.dat','r');
dTdH = fscanf(fileID, '%f', [1,Inf]);
fileID = fopen('build/dTdC.dat','r');
dTdC = fscanf(fileID, '%f', [1,Inf]);

seed = 101;
fclose(fileID);

CD = reshape(CD, seed, seed);
HD = reshape(HD, seed, seed);
TD = reshape(TD, seed, seed);
Vf = reshape(Vf, seed, seed);

dTdH = reshape(dTdH, seed, seed);
dTdC = reshape(dTdC, seed, seed);

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

figure
g = surf(CD, HD, dTdH);
set(g,'linestyle','none','facecolor','interp');
light("Style","local","Position",[0 0 10]);
title("Composition-Enthalpy-dTdH");
xlabel("Composition");
ylabel("Enthalpy");
zlabel("dTdH");

figure
g = surf(CD, HD, dTdC);
set(g,'linestyle','none','facecolor','interp');
light("Style","local","Position",[0 0 10]);
title("Composition-Enthalpy-dTdC");
xlabel("Composition");
ylabel("Enthalpy");
zlabel("dTdC");


% ============================================================
