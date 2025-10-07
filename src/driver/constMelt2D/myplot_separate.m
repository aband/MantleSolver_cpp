function [] = myplot_separate(M, N, folder, k)

filename = strcat(folder, '/gridCellX.dat');
fileID = fopen(filename, 'r');
pX = fscanf(fileID, '%f', [1, Inf]);

filename = strcat(folder, '/gridCellY.dat');
fileID = fopen(filename, 'r');
pY = fscanf(fileID, '%f', [1, Inf]);

pX = reshape(pX, M, N);
pY = reshape(pY, M, N);

%v = VideoWriter('videoall.avi','Motion JPEG AVI');
%open(v);

%fullname = strcat(folder, '/porosity');
%fullname = strcat(fullname, '*.dat');
%fstruct1 = dir(fullname);
%fcell1 = struct2cell(fstruct1);

%loops = numel(fstruct1)

%h = figure;

cut = 80 - 80*1.8/2 + 1;

time = k*2*10;
% Porosity
figure
filename = strcat(folder, '/porosity');
filename = strcat(filename, string(k));
filename = strcat(filename, '.dat');
fileID = fopen(filename, 'r');
porosity = fscanf(fileID, '%f ', [1, Inf]);
porosity = reshape(porosity, M, N);

%subplot(1,4,1)
plot(porosity(2,:), pY(2,:), 'LineWidth',3);
xlim([0,0.15])
mytitle = strcat('Porosity, t = ', string(time));
t = title(mytitle)
ylabel('Depth');
set(gcf, 'Position',[50 50 400 1200]);
fontsize(16,"points")

figure
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

plot(-darcypressure(2,:), pY(2,:), 'LineWidth',3);
hold on
plot(stokespressure(2,:), pY(2,:),'--','LineWidth',3);
hold off
mytitle = strcat('Pressure, t = ', string(time));
t = title(mytitle)
ylabel('Depth');
legend('Darcy','Stokes');
set(gcf, 'Position',[50 50 400 1200]);
fontsize(16,"points")

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

figure
plot(unscaleddarcyy(2,:), pY(2,:),'LineWidth',3);
hold on
plot(stokesy(2,:), pY(2,:),'--','LineWidth',3);
hold off
mytitle = strcat('Velocity, t = ', string(time));
t = title(mytitle)
xlim([-5e-3, 5e-3])
ylabel('Depth')
legend('Darcy','Stokes','Location','southeast')
set(gcf, 'Position',[50 50 400 1200]);
fontsize(16,"points")
