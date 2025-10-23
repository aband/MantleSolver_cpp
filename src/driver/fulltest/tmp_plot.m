grid = [1.0/8, 1.0/16, 1.0/32, 1.0/64]
grid = [8, 16, 32, 64];

figure
loglog(grid, [.4502 , .1052, .02513, .006130],'linewidth',2);
hold on
loglog(grid, [1.039, .4463, .2093, .1019],'linewidth',2);
loglog(grid, [2.632, .5845, .1374 .03330],'linewidth',2);

legend('u', 'du', 'p')
title('Convergence on Rectangular Grids')
xlabel('grid')
ylabel('error')

figure
loglog(grid, [.4603 .1074 .02567 .006266],'linewidth',2);
hold on
loglog(grid, [1.070 .4670 .2207 .1078],'linewidth',2);
loglog(grid, [2.685 .5981 .1525 .04580],'linewidth',2);

legend('u', 'du', 'p')
title('Convergence on Distorted Grids')
xlabel('grid')
ylabel('error')


