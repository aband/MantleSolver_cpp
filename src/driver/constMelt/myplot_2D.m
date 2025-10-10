function [] = myplot_all(M, N, folder)

filename = strcat(folder, '/gridCellX.dat');
fileID = fopen(filename, 'r');
pX = fscanf(fileID, '%f', [1, Inf]);

filename = strcat(folder, '/gridCellY.dat');
fileID = fopen(filename, 'r');
pY = fscanf(fileID, '%f', [1, Inf]);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

v = VideoWriter('videoall.avi','Motion JPEG AVI');
open(v);

fullname = strcat(folder, '/porosity');
fullname = strcat(fullname, '*.dat');
fstruct1 = dir(fullname);
fcell1 = struct2cell(fstruct1);

loops = numel(fstruct1)

h = figure;

set(gcf, 'Position',[50 50 1800 700]);

cut = 80 - 80*1.8/2 + 1;

for k=1:loops

filename = strcat(folder, '/porosity');
filename = strcat(filename, string(k));
filename = strcat(filename, '.dat');
fileID = fopen(filename, 'r');
porosity = fscanf(fileID, '%f ', [1, Inf]);
porosity = reshape(porosity, M, N);

subplot(1,3,1) 
%plot(porosity(2,:), pY(2,:), 'LineWidth',3);
surf(pX, pY, porosity)
colormap turbo
shading interp
%xlim([0,0.1])
title('porosity')
ylabel('Depth');

filename = strcat('build/stokesVx',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
stokesx = fscanf(fileID, '%f', [1,Inf]);

filename = strcat('build/stokesVy',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
stokesy = fscanf(fileID, '%f', [1,Inf]);

filename = strcat('build/darcyVx',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
darcyx = fscanf(fileID, '%f', [1,Inf]);

filename = strcat('build/darcyVy',string(k));
filename = strcat(filename,'.dat');
fileID = fopen(filename, 'r');
darcyy = fscanf(fileID, '%f', [1,Inf]);

%filename = strcat('build/porosity',string(k));
%filename = strcat(filename,'.dat');
%fileID = fopen(filename, 'r');
%poro   = fscanf(fileID, '%f', [1,Inf]);
%
stokesx = reshape(stokesx, M, N);
stokesy = reshape(stokesy, M, N);
darcyx = reshape(darcyx, M, N);
darcyy = reshape(darcyy, M, N);
%poro   = reshape(poro,M,N);

subplot(1,3,2)
quiver(pX, pY, stokesx, stokesy);
title(["Stokes Velocity",num2str(k)])

unscaleddarcyy = darcyy.*porosity;

%subplot(2,3,3)
subplot(1,3,3)
quiver(pX, pY, darcyx.*porosity, unscaleddarcyy);
title(["Darcy Velocity",num2str(k)])

%{
filename = strcat(folder, '/qf');
filename = strcat(filename, string(k));
filename = strcat(filename, '.dat');
fileID = fopen(filename, 'r');
darcypressure = fscanf(fileID, '%f ', [1, Inf]);
darcypressure = reshape(darcypressure, M, N);

filename = strcat(folder, '/qs');
filename = strcat(filename, string(k));
filename = strcat(filename, '.dat');
fileID = fopen(filename, 'r');
stokespressure = fscanf(fileID, '%f ', [1, Inf]);
stokespressure = reshape(stokespressure, M, N);

subplot(1,4,3)
plot(-darcypressure(2,:), pY(2,:), 'LineWidth',3);
hold on
plot(stokespressure(2,:), pY(2,:),'--','LineWidth',3);
hold off
title('Pressure')
ylabel('Depth');
legend('Darcy Pressure','Stokes Pressure');

subplot(1,4,4)
plot(stokespressure(2,cut:end)+darcypressure(2,cut:end), pY(2,cut:end), 'LineWidth', 3);
title('Effective Pressure')
ylabel('Depth')

filename = strcat(folder, '/stokesVx',string(k));
filename = strcat(filename, '.dat');
fileID = fopen(filename, 'r');
stokesx = fscanf(fileID, '%f', [1, Inf]);

filename = strcat(folder, '/stokesVy',string(k));
filename = strcat(filename, '.dat');
fileID = fopen(filename, 'r');
stokesy = fscanf(fileID, '%f', [1, Inf]);

stokesx = reshape(stokesx, M, N);
stokesy = reshape(stokesy, M, N);

filename = strcat(folder, '/darcyVx',string(k));
filename = strcat(filename, '.dat');
fileID = fopen(filename, 'r');
darcyx = fscanf(fileID, '%f', [1, Inf]);

filename = strcat(folder, '/darcyVy',string(k));
filename = strcat(filename, '.dat');
fileID = fopen(filename, 'r');
darcyy = fscanf(fileID, '%f', [1, Inf]);

darcyx = reshape(darcyx, M, N);
darcyy = reshape(darcyy, M, N);

unscaleddarcyy = darcyy.*porosity;

subplot(1,4,2)
plot(stokesy(2,:), pY(2,:),'LineWidth',3);
hold on
plot(unscaleddarcyy(2,:), pY(2,:),'--','LineWidth',3);
hold off
title('Velocity')
%xlim([-5e-3, 5e-3])
ylabel('Depth')
legend('Stokes Vel','Darcy Vel','Location','southeast')
%}
time = k*10*1

mytitle = strcat('Time =  ', string(time));

sgtitle(mytitle);

pause
G = getframe(gcf);

writeVideo(v,G);

end

close(v);
