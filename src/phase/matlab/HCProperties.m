% file: HCDiagrams.m
% date: 24 Jul 2019
% author: Marc A. Hesse

% Description: Plot all inverted quanties as functions (tD and phi's) as
% function of CD and HD

close all, clear, clc

%Phase.L  = 3.34e5;
%Phase.Xe = 0.8859;
%Phase.T1 = 273.15; % melting temperature of phase 1
%Phase.Te = 252.05; % eutectic temperature

Phase.L  = 4e5;
Phase.Xe = 0.7;
Phase.T1 = 1350;
Phase.Te = 1227;

[rho,cp,kappa,hD,HD,Phase] = physical_properties_eutectic(Phase);
[HD_of_TX,CD_of_TX,phi_of_TX,f_of_TX,Phase] = setup_phase_behavior(rho,cp,hD,Phase);

CDmax = max(Phase.invHC.CD3l,Phase.invHC.CD3s)
Phase.invHC.CD3l
Phase.invHC.CD3s

HDmax = 15;

Nc = 10; Nh = 10;
cc = linspace(0,CDmax,Nc);
hh = linspace(-1,HDmax,Nh);
[CD,HD] = meshgrid(cc,hh);

[TD,Phi,reg] = eval_phase_behavior(HD,CD,rho,cp,hD,Phase);

[pTD,pPhi,preg] = eval_phase_behavior(6,0.2,rho,cp,hD,Phase)

TD = reshape(TD,Nh,Nc);
Phi_ice = reshape(Phi(:,1),Nh,Nc);
Phi_sal = reshape(Phi(:,2),Nh,Nc);
Phi_bri = reshape(Phi(:,3),Nh,Nc);

scsz = get(0,'ScreenSize');
scw = scsz(3); sch = scsz(4);
aspect = 1/0.8;

% figure size
figW = scw/2;
gapL = 0.06*figW;
gapW = 0.06*figW;
gapR = 0.02*figW;
gapB = 0.05*figW;
gapT = 0.02*figW;
subW = (figW - (gapL+gapR+3*gapW))/4;
subH = subW*aspect;
figH = gapB+subH+gapT;

% subplot locations
sub1x = gapL/figW;
sub1y = (gapB)/figH;

sub2x = sub1x + (subW+gapW)/figW;
sub2y = sub1y;

sub3x = sub2x + (subW+gapW)/figW;
sub3y = sub1y;

sub4x = sub3x + (subW+gapW)/figW;
sub4y = sub1y;

figure('Position',[0, 0, figW, figH])
set(gca,'FontSize',70)
subplot('position',[sub1x sub1y subW/figW subH/figH])
contourf(CD,HD,TD,30), hold on, colorbar('location','east')
plot_HC_diagram(Phase,CDmax,-1,HDmax)
title('T_D')

subplot('position',[sub2x sub2y subW/figW subH/figH])
contourf(CD,HD,Phi_ice,30), hold on, colorbar('location','east')
plot_HC_diagram(Phase,CDmax,-1,HDmax)
caxis([0 1])

title('\phi_{ice}')

subplot('position',[sub3x sub3y subW/figW subH/figH])
contourf(CD,HD,Phi_sal,30), hold on, colorbar('location','east')
plot_HC_diagram(Phase,CDmax,-1,HDmax)
caxis([0 1])

title('\phi_{sal}')

subplot('position',[sub4x sub4y subW/figW subH/figH])
contourf(CD,HD,Phi_bri,30), hold on, colorbar('location','east')
caxis([0 1])
plot_HC_diagram(Phase,CDmax,-1,HDmax)
title('\phi_{bri}')
