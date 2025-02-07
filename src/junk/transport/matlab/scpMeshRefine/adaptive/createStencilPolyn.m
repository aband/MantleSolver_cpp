function [basepolyncoeff32All, basepolyncoeff43All] = createStencilPolyn(grid)

N = length(grid)-1;

global stencil32
global linWgt32  
global stencil32L
global linWgt32L 
global stencil32R 
global linWgt32R  
global stencil43
global linWgt43  
global stencil43L 
global linWgt43L  
global stencil43LL
global linWgt43LL  
global stencil43R
global linWgt43R   
global nStencils32
global nStencils43 
global maxR32 
global maxR43 

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
