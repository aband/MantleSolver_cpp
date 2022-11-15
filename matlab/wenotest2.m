function wenotest2(totLength, totTime, iMax, nMax, plotUBar, runType, fluxFcn)

stencil = [ [-2,2]; [-2,0]; [-1,1]; [0,2]; [0,0] ];
linWgt  = [ 1/2 1/8 1/8 1/8 1/8];
%stencil = [ [-2,0]; [-1,1]; [0,2] ];
%linWgt  = [ 1/3 1/3 1/3];
%stencil = [ [-2,2]];
%linWgt  = [ 1 ];


nStencils = size(stencil,1);
if nStencils ~= length(linWgt)
    warning('Wrong number of linear weights.');
    return;
end

for s=1:nStencils
    stencilDegree = stencil(s,2) - stencil(s,1) + 1;
    switch stencilDegree
        case 1
        %case 2
        case 3
        case 5
        %case 7
      otherwise
        warning('Unsupported stencil length.');
        return;
    end
    if linWgt < 0
        warning('Negative linear weight.');
        return;
    end
end

if sum(linWgt) == 0
    warning('All zero linear weights.');
    return;
end
linWgt = linWgt / sum(linWgt);

% FLUX_FUNCTION

switch fluxFcn
  case 'linear'
    f = @(u) u;
    df = @(u) 1;
    LF = 1;
  case 'burgers'
    f = @(u) u.^2/2;
    df = @(u) u;
    LF = 1;
  case 'BL' %Buckley-Leverett
    f = @(u) u^2 / ( u^2 + (1-u)^2 );
    df = @(u) 2*u*(1-u) / ( u^2 + (1-u)^2 )^2;
    LF = 2;
  otherwise
    warning('Unexpected flux function type. Exiting.');
    return;
end

% RUN_TYPE: IC

trueUKnown = 0;

switch runType
  %%Step Down
  case 'stepDown'
    icVal = @(x) x<=0.5*totLength;
    switch fluxFcn
      case 'linear'
        trueUKnown = 1;
        trueU = @(x,t) mod(x-t,totLength)<=0.5*totLength;
      case 'burgers'
        trueUKnown = 1;
        trueU = @(x,t) mod(x-t/2,totLength)<=0.5*totLength;
      otherwise
    end

  %%Step Up
  case 'stepUp'
    icVal  = @(x) x>=0.5*totLength;
    switch fluxFcn
      case 'linear'
        trueUKnown = 1;
        trueU = @(x,t) mod(x-t,totLength)>=0.5*totLength;
      case 'burgers'
        trueUKnown = 1;
        trueU = @(x,t) max(0,min((x-0.5*totLength)/t,1));
      otherwise
        trueUKnown = 0;
    end
 
  %%Sine wave
  case 'sine'
    icVal  = @(x) 0.5*( 1 + sin(2*pi*(x)) ) ;
    switch fluxFcn
      case 'linear'
        trueUKnown = 1;
        trueU = @(x,t) 0.5*( 1 + sin(2*pi*(x-t)) ) ;
      otherwise
        trueUKnown = 0;
    end
    
  %%Saw tooth
  case 'saw'
    power = 1;
    icVal  = @(x) x^power;
    switch fluxFcn
      case 'linear'
        trueUKnown = 1;
        trueU = @(x,t) mod((x-t),totLength)^power;
      otherwise
        trueUKnown = 0;
    end
  otherwise
    warning('Unexpected run (IC) type. Exiting.');
    return;
end

% INITIALIZE

hatF = @(a,b) 0.5*( f(a) + f(b) - LF*(b-a) );

gaussPt = [-sqrt(3/5) 0 sqrt(3/5)];
gaussWt = [5/18 8/18 5/18]; %half value!

format compact;

pauseOutput = 0;
plotTime = 6;

if totLength <= 0
    totLength = 1;
end
if totTime <= 0
    totTime = 1;
end
if iMax < 5
    iMax = 5;
end
if nMax < 1
    nMax = 1;
end

uRK1 = zeros(iMax,1);
uRK2 = zeros(iMax,1);
uBar = zeros(iMax,1);
if trueUKnown
    uTrue = zeros(iMax,1);
end

dt = totTime/nMax;
dx = totLength/iMax;

lambda = dt/dx;
if lambda*LF > 1
    warning('CFL constraint violated.');
end

for i=1:1:iMax
    xL = (i-1)*dx;
    xR = i*dx;
    gPt = xL + (xR-xL)*(gaussPt+1)/2;
    uBar(i) = gaussWt*[icVal(gPt(1)); icVal(gPt(2)); icVal(gPt(3))];
end

% COMPUTE_SOLUTION

for n=1:1:nMax
    
    time = n*dt;
 
    % RK Time Step Stage 1
    
    [uL,uR] = reconstructionLR(iMax,dx,uBar,stencil,linWgt);
    for i=1:1:iMax
        iL = 1 + mod(i-2,iMax);
        iR = 1 + mod(i  ,iMax);
        uRK1(i) = uBar(i) - lambda*( hatF(uR(i),uL(iR)) - hatF(uR(iL),uL(i)) );
    end
    
    % RK Time Step Stage 2
    
    [uL,uR] = reconstructionLR(iMax,dx,uRK1,stencil,linWgt);
    for i=1:1:iMax
        iL = 1 + mod(i-2,iMax);
        iR = 1 + mod(i  ,iMax);
        uRK2(i) = (3*uBar(i) + uRK1(i))/4 - (1/4)*lambda*( hatF(uR(i),uL(iR)) - hatF(uR(iL),uL(i)) );
    end
    
    % RK Time Final Stage
    
    [uL,uR] = reconstructionLR(iMax,dx,uRK2,stencil,linWgt);
    for i=1:1:iMax
        iL = 1 + mod(i-2,iMax);
        iR = 1 + mod(i  ,iMax);
        uBar(i) = (uBar(i) + 2*uRK2(i))/3 - (2/3)*lambda*( hatF(uR(i),uL(iR)) - hatF(uR(iL),uL(i)) );
    end
    
    % Mass consevation
    %totMass = sum(uBar)*dx;
    %fprintf('Total Mass: %g\n',totMass);

    % Plot Solution
    
    if plotUBar
        clf;

        x = [0 dx/2:dx:totLength-dx/2 totLength];

        plot(x,ones(iMax+2,1),'g','LineWidth',2);
        hold on;
        plot(x,zeros(iMax+2,1),'g','LineWidth',2);
            
        if trueUKnown
            for i=1:1:iMax
                xL = (i-1)*dx;
                xR = i*dx;
                gPt = ( xL + (xR-xL)*(gaussPt+1)/2 );
                uTrue(i) = gaussWt*[trueU(gPt(1),time); trueU(gPt(2),time); trueU(gPt(3),time)];
            end
            uAvg = 0.5*(uTrue(1) + uTrue(iMax));
            plot(x,[uAvg; uTrue; uAvg],'k','LineWidth',2);
        end

        uAvg = 0.5*(uBar(1) + uBar(iMax));
        plot(x,[uAvg; uBar; uAvg],'b','LineWidth',2);

        ax = gca;
        ax.FontSize = 20;
        ax.LineWidth = 2;
        ax.XLim = [0 totLength];
        ax.YLim = [-0.2 1.2];
            
        if pauseOutput
            pause;
        else
            % plot in plotTime seconds
            pause(min(1,plotTime/nMax));
        end
        hold off;
    end
end

% COMPUTE_ERRORS

if trueUKnown
    error = sum(abs(uBar - uTrue))*dx;

    fprintf('Computed L1 Error: %5.3e\n', error);
end

end


function [uL,uR] = reconstructionLR(iMax,dx,uBar,stencil,linWgt)
    eps0 = 1;
    
    nStencils = size(stencil,1);

    % Nonlinear weights

    hatWgt = linWgt;
    % -----------------------------------

    uL = zeros(iMax,1);
    uR = zeros(iMax,1);
    
    maxR = max(stencil(:,2) - stencil(:,1) + 1);
    basePolynCoeff = zeros(nStencils,maxR,maxR);

    for s=1:nStencils
        r = stencil(s,2) - stencil(s,1) + 1;
        basePolynCoeff(s,1:r,1:r) = polyn(stencil(s,:));
    end
    
    for i=1:1:iMax

        iLLL = 1 + mod(i-4,iMax);
        iLL  = 1 + mod(i-3,iMax);
        iL   = 1 + mod(i-2,iMax);
        iR   = 1 + mod(i  ,iMax);
        iRR  = 1 + mod(i+1,iMax);
        iRRR = 1 + mod(i+2,iMax);

        fullStencil = [iLLL,iLL,iL,i,iR,iRR,iRRR];

        for s=1:nStencils

            left = stencil(s,1)+4;
            right = stencil(s,2)+4;

            r = right - left + 1;

            currentStencil = fullStencil(left:right);

            uBarStencil = uBar(currentStencil);

            % Create smoothness indicator

%{
 {            dis = (left-4):(right-4);
 {            uBar_diff = uBarStencil - uBarStencil(dis==0);
 {            uBar_diff = uBar_diff(dis~=0);
 {            dis = dis(dis~=0);
 {
 {            uBar_diff = (uBar_diff.*uBar_diff)'./(dis.*dis);
 {
 {            sigma = 1/(size(uBar_diff,2))*sum(uBar_diff);
 %}

            sigma  = computeSigma(uBarStencil,stencil(s,:));

            r = stencil(s,2) - stencil(s,1) + 1;
            eta = floor(r/2)+1;

            hatWgt(s) = linWgt(s) / ( sigma^eta + eps0*(dx/100)^r);

       end 

       nonlinWgt = hatWgt / sum(hatWgt)

       for s=1:nStencils

            left = stencil(s,1)+4;
            right = stencil(s,2)+4;

            currentStencil = fullStencil(left:right);

            uBarStencil = uBar(currentStencil);

            % Create coefficients of polynomial
            r = right - left + 1;

            uL(i) = uL(i) + nonlinWgt(s)*weno_reconst(-0.5,basePolynCoeff,s,uBarStencil,r); 
            uR(i) = uR(i) + nonlinWgt(s)*weno_reconst( 0.5,basePolynCoeff,s,uBarStencil,r); 
       end

    end
end

function [ru] = weno_reconst(p,coeff,s,uBarStencil,r)

   % p is relative coordinate (x-x0)/dx

   ru = 0.0;

   for k=1:r
       for i = 1:r
           ru = ru + uBarStencil(k) * coeff(s,k,i)*p.^(i-1);
       end
   end

end

function [sol] = polyn(stencil)

r = stencil(2) - stencil(1) + 1;

M = zeros(r,r);

sol = zeros(r,r);

    for k = 1:r
       for j = 1:r
           xleft  = stencil(1)+(j-1)-0.5;
           xright = xleft+1;

           for i=1:r
               M(j,i) = integral(@(x) x.^(i-1), xleft, xright);
           end
       end

       B = zeros(r,1);
       B(k) = 1;
       sol(k,:) = M\B;
    end

end

function [Sigma] = computeSigma(uBar, stencil)

    Sigma = 0;

    target = 1-stencil(1);

    if stencil(2)> stencil(1)
        for i = stencil(1):stencil(2)
            if i~= 0
                I = i-stencil(1) + 1;
                Sigma = Sigma + (uBar(I) - uBar(target))^2;
                %Sigma = Sigma + ( ( uBar(I) - uBar(target) ) / i )^2;
            end
        end
        %Sigma = Sigma /(stencil(2) - stencil(1));
        %midIndex = floor((stencil(1) + stencil(2))/2) - stencil(1) + 1;
        %Sigma = Sigma + ( uBar(1) + uBar(end) - 2*uBar(midIndex) )^2;
 
    else
        Sigma = 0;
    end

end

function [Sigma] = computeSigma2(uBar, stencil)

    Sigma = 0;

    target = 1-stencil(1);

    if stencil(2)> stencil(1)
        for i = stencil(1):stencil(2)
            %if i~= 0
						  for j = stencil(1):stencil(2)
                I = i-stencil(1) + 1;
					 J = j-stencil(1) + 1;
                %Sigma = Sigma + (uBar(I) - uBar(target)^2);
                %Sigma = Sigma + ( ( uBar(I) - uBar(target) ) / i )^2;
					 if J ~= I
				    Sigma = Sigma + ( (uBar(I) - uBar(J) ) / (J-I) )^2;
					 end
						  end
            %end
        end
		  Sigma = Sigma / ((stencil(2)-stencil(1))*(stencil(2)-stencil(1)+1)/2);
        %Sigma = Sigma /(stencil(2) - stencil(1));
        %midIndex = floor((stencil(1) + stencil(2))/2) - stencil(1) + 1;
        %Sigma = Sigma + ( uBar(1) + uBar(end) - 2*uBar(midIndex) )^2;
 
    else
        Sigma = 0;
    end

end


