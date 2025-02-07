function [grid] = generateMesh(a,N)

b = 2/(N*(N-1))*(1-a*N);

multi = 0:(N-1);

interval = multi*b+a;

grid = [0];

for i = 1:N

grid = [grid,sum(interval(1:i))];

end

grid = 1-grid;

grid = -grid;

grid = [grid(1:end-1),-flip(grid)];

end
