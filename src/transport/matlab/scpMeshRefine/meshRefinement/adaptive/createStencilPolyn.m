function [basepolyncoeff32All, basepolyncoeff43All] = createStencilPolyn(grid)

N = length(grid)-1;

stencil32 = [[-1,1]; [-1,0]; [0,1]];
linWgt32  = [3,1,1];

stencil32L = [[0,2];[0,1];[0,0]];
linWgt32L  = [3,2,1];

stencil32R = [[-2,0];[-1,0];[0,0]];
linWgt32R  = [3,2,1];

% Weno (4,3) reconstruction for diffusive flux
stencil43 = [[-2,1];[-2,0];[-1,1];[0,0]];
linWgt43  = [4,1,1,1]; 

stencil43L = [[0,3];[0,2];[0,1];[0,0]];
linWgt43L  = [4,1,1,1]; 

stencil43LL = [[-1,2];[-1,1];[0,2];[0,0]];
linWgt43LL  = [4,1,1,1]; 

stencil43R = [[-3,0];[-2,0];[-1,0];[0,0]];
linWgt43R  = [4,1,1,1]; 

% Create stencil polynomials with given grid information

nStencils32 = size(stencil32,1);
nStencils43 = size(stencil43,1); 

maxR32 = max(stencil32(:,2)-stencil32(:,1)+1);
maxR43 = max(stencil43(:,2)-stencil43(:,1)+1);

basepolyncoeff32All = zeros(N,nStencils32,maxR32,maxR32);
basepolyncoeff43All = zeros(N+1,nStencils43,maxR43,maxR43);

basepolyncoeff32L = basePolynRefine(1,stencil32L,-0.5,grid);
basepolyncoeff43L = basePolynRefine(1,stencil43L, 0.0,grid);
basepolyncoeff43LL = basePolynRefine(2,stencil43LL, 0.0,grid);

basepolyncoeff32R = basePolynRefine(N,stencil32R,-0.5,grid);
basepolyncoeff43R = basePolynRefine(N,stencil43R, 0.0,grid);
basepolyncoeff43RR = basePolynRefine(N,stencil43R, -1.0,grid);

basepolyncoeff32All(1,:,:,:) = basepolyncoeff32L;
basepolyncoeff32All(N,:,:,:) = basepolyncoeff32R;

for s = 2:N-1
    basepolyncoeff32All(s,:,:,:) = basePolynRefine(s,stencil32,-0.5,grid);
end

basepolyncoeff43All(1,:,:,:)   = basepolyncoeff43L;
basepolyncoeff43All(2,:,:,:)   = basepolyncoeff43LL;
basepolyncoeff43All(N,:,:,:)   = basepolyncoeff43RR;
basepolyncoeff43All(N+1,:,:,:) = basepolyncoeff43R;

for s = 3:N-1
    basepolyncoeff43All(s,:,:,:) = basePolynRefine(s,stencil43, 0.0,grid);
end

end
