function [] = order_plot(M, N, mark)

fileID = fopen('build/eff_order1.dat','r');
order = fscanf(fileID, '%f ', [1, Inf]);

order = reshape(order, M, N)

px = 1:M;
py = 1:N;

[pX, pY] = meshgrid(px, py);

order = [order, ones(N,1);
         ones(1,M), 1];

figure
s = pcolor(order'-1); 
%colorbar
clim([1.8, 4.2])
colormap(turbo)
axis equal
title('Effective Order')

fig = figure 
ax = axes;
colormap(turbo)
caxis([1.8,4.2])
c = colorbar(ax)
ax.Visible = 'off'

exportgraphics(fig, 'colorbar.png','resolution',300)
