% Implement
%    u_t + f(u)_x = 0
% using periodic BCs, the classic WENO(5,3) FV scheme,
% and three stage explicit SSP-RK3 time stepping:
%   v_0 = u^n
%   v_1 = u^n + dt F(u^n)
%   v_2 = (3/4) u^n + (1/4) v^1 + (1/4) dt F(v^1)
%   u^{n+1} = (1/3) u^n + (2/3) v^2 + (2/3) dt F(v^2)
% requiring CFL <= 1.
%
% Call as: weno53(totLength, totTime, iMax, nMax, plotUBar, runType, fluxFcn);
function weno53(totLength, totTime, iMax, nMax, plotUBar, runType, fluxFcn)

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
    icVal  = @(x) 0.5*( 1 + sin(2*pi*(x)) );
    switch fluxFcn
      case 'linear'
        trueUKnown = 1;
        trueU = @(x,t) 0.5*( 1 + sin(2*pi*(x-t)) );
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
    
    [uL,uR] = reconstructionLR(iMax,dx,uBar);
    for i=1:1:iMax
        iL = 1 + mod(i-2,iMax);
        iR = 1 + mod(i  ,iMax);
        uRK1(i) = uBar(i) - lambda*( hatF(uR(i),uL(iR)) - hatF(uR(iL),uL(i)) );
    end
    
    % RK Time Step Stage 2
    
    [uL,uR] = reconstructionLR(iMax,dx,uRK1);
    for i=1:1:iMax
        iL = 1 + mod(i-2,iMax);
        iR = 1 + mod(i  ,iMax);
        uRK2(i) = (3*uBar(i) + uRK1(i))/4 - (1/4)*lambda*( hatF(uR(i),uL(iR)) - hatF(uR(iL),uL(i)) );
    end
    
    % RK Time Final Stage
    
    [uL,uR] = reconstructionLR(iMax,dx,uRK2);
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


function [uL,uR] = reconstructionLR(iMax,dx,uBar)
    uL = zeros(iMax,1);
    uR = zeros(iMax,1);

    for i=1:1:iMax
        iLL = 1 + mod(i-3,iMax);
        iL  = 1 + mod(i-2,iMax);
        iR  = 1 + mod(i  ,iMax);
        iRR = 1 + mod(i+1,iMax);

        % Values of the polynomials and smoothness indicators
        [val3L,val3R] = polyn3(uBar,iLL,iL,i,iR,iRR);
        sigma3 = smoothness3(uBar,iLL,iL,i,iR,iRR);

        [val5L,val5R] = polyn5(uBar,iLL,iL,i,iR,iRR);

        % Linear weights
        alphaL = [0.3 0.6 0.1];
        alphaR = [0.1 0.6 0.3];
        
        % Nonlinear weights
        alphaL = alphaL./((sigma3 + (1e-4)*dx^2).^2);
        alphaR = alphaR./((sigma3 + (1e-4)*dx^2).^2);
        alphaL = alphaL / ( alphaL(1) + alphaL(2) + alphaL(3) );
        alphaR = alphaR / ( alphaR(1) + alphaR(2) + alphaR(3) );
    
        % Values of the reconstruction
        %uL(i) = alphaL*val3L';
        %uR(i) = alphaR*val3R';

        uL(i) = val5L';
        uR(i) = val5R';
    end
end


function [val3M,val3P] = polyn3(u,iLL,iL,i,iR,iRR)
  val3M(1) = (-1/6)*u(iLL) + ( 5/6)*u(iL) + ( 1/3)*u(i);
  val3P(1) = ( 1/3)*u(iLL) + (-7/6)*u(iL) + (11/6)*u(i);

  val3M(2) = ( 1/3)*u(iL) + (5/6)*u(i) + (-1/6)*u(iR);
  val3P(2) = (-1/6)*u(iL) + (5/6)*u(i) + ( 1/3)*u(iR);
  
  val3M(3) = (11/6)*u(i) + (-7/6)*u(iR) + ( 1/3)*u(iRR);
  val3P(3) = ( 1/3)*u(i) + ( 5/6)*u(iR) + (-1/6)*u(iRR);
end


function [val5M,val5P] = polyn5(u,iLL,iL,i,iR,iRR)
  val5M = (-1/20)*u(iLL) + (  9/20)*u(iL) + (47/60)*u(i) + (-13/60)*u(iR) + ( 1/30)*u(iRR);
  val5P = ( 1/30)*u(iLL) + (-13/60)*u(iL) + (47/60)*u(i) + (  9/20)*u(iR) + (-1/20)*u(iRR);
end


function [sigma3] = smoothness3(u,iLL,iL,i,iR,iRR)
  sigma3(1) = (13/12)*( u(iLL) - 2*u(iL) + u(i  ) )^2 + (1/4)*( u(iLL) - 4*u(iL) + 3*u(i  ) )^2;
  sigma3(2) = (13/12)*( u(iL ) - 2*u(i ) + u(iR ) )^2 + (1/4)*( u(iL) - u(iR) )^2;
  sigma3(3) = (13/12)*( u(i  ) - 2*u(iR) + u(iRR) )^2 + (1/4)*( 3*u(i) - 4*u(iR) + u(iRR) )^2;
end
