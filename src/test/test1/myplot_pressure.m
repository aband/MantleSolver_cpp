function [] = myplot_pressure(M, N)

fileID = fopen('build/gridCellX.dat','r');
pX = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen('build/gridCellY.dat','r');
pY = fscanf(fileID, '%f', [1,Inf]);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

fileID = fopen('build/referencep.dat','r');
referencep = fscanf(fileID, '%f', [1,Inf]);

fstruct1 = dir('build/*stokesq*.dat');
fstruct2 = dir('build/*darcyq*.dat');
fstruct3 = dir('build/*stokesp*.dat');
fstruct4 = dir('build/*darcyp*.dat');

fstruct1 = rmfield(fstruct1,'folder');
fstruct1 = rmfield(fstruct1,'date');
fstruct1 = rmfield(fstruct1,'bytes');
fstruct1 = rmfield(fstruct1,'isdir');
fstruct1 = rmfield(fstruct1,'datenum');

fstruct2 = rmfield(fstruct2,'folder');
fstruct2 = rmfield(fstruct2,'date');
fstruct2 = rmfield(fstruct2,'bytes');
fstruct2 = rmfield(fstruct2,'isdir');
fstruct2 = rmfield(fstruct2,'datenum');

fstruct3 = rmfield(fstruct3,'folder');
fstruct3 = rmfield(fstruct3,'date');
fstruct3 = rmfield(fstruct3,'bytes');
fstruct3 = rmfield(fstruct3,'isdir');
fstruct3 = rmfield(fstruct3,'datenum');

fstruct4 = rmfield(fstruct4,'folder');
fstruct4 = rmfield(fstruct4,'date');
fstruct4 = rmfield(fstruct4,'bytes');
fstruct4 = rmfield(fstruct4,'isdir');
fstruct4 = rmfield(fstruct4,'datenum');

fcell1 = struct2cell(fstruct1);
fcell2 = struct2cell(fstruct2);
fcell3 = struct2cell(fstruct3);
fcell4 = struct2cell(fstruct4);

for k=1:numel(fstruct1)

fileID = fopen(strcat('build/',fcell1{k}), 'r');
stokesq = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen(strcat('build/',fcell2{k}), 'r');
darcyq = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen(strcat('build/',fcell3{k}), 'r');
stokesp = fscanf(fileID, '%f', [1,Inf]);

fileID = fopen(strcat('build/',fcell4{k}), 'r');
darcyp = fscanf(fileID, '%f', [1,Inf]);

stokesq = reshape(stokesq, M, N);
darcyq = reshape(darcyq, M, N);
stokesp = reshape(stokesp, M, N);
darcyp = reshape(darcyp, M, N);

figure
subplot(2,2,1)
surf(pX, pY, stokesq);
title(["Stokes Pressure Potential",num2str(k)])

subplot(2,2,2)
surf(pX, pY, darcyq);
title(["Darcy Pressure Potential",num2str(k)])

subplot(2,2,3)
surf(pX, pY, stokesp);
title(["Stokes Pressure",num2str(k)])

subplot(2,2,4)
surf(pX, pY, darcyp);
title(["Darcy Pressure",num2str(k)])

figure
subplot(2,2,1)
plot(stokesq(1,:),pY(2,:));
title(["Stokes Pressure Potential",num2str(k)])
ylabel("Depth");

subplot(2,2,2)
plot(darcyq(1,:),pY(2,:));
title(["Darcy Pressure Potential",num2str(k)])
ylabel("Depth")

subplot(2,2,3)
plot(stokesp(1,:), pY(2,:));
title(["Stokes Pressure",num2str(k)])
ylabel("Depth")

subplot(2,2,4)
plot(darcyp(1,:),pY(2,:));
title(["Darcy Pressure",num2str(k)])
ylabel("Depth")

end

fclose(fileID);
