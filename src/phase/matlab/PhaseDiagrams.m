% file test_quadraticCH.m
% file: PhaseDiagrams.m
% author: Marc Hesse
% date: 19 Jul 2019
% Description: Make figure showing all the different ways of plotting the
% phase diagram for supplementary materials.

close all, clear, clc
%Phase.L  = 3.34e5;
%Phase.Xe = 0.8859;
%Phase.T1 = 273.15;
%Phase.Te = 252.05;

Phase.L  = 4e5;
Phase.Xe = 0.7;
Phase.T1 = 1350;
Phase.Te = 1227;

% [rho,cp,kappa,hD,Phase] = physical_properties_eutectic(Phase);
[rho,cp,kappa,hD,HD_from_TPhi,Phase] = physical_properties_eutectic(Phase);

Phase

rho

cp

kappa

hD

HD_from_TPhi

[HD_of_TX,  CD_of_TX,  phi_of_TX,f_of_TX,Phase] = setup_phase_behavior(rho,cp,hD,Phase);

HD_of_TX

CD_of_TX

phi_of_TX

Phase

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

subplot('position',[sub1x sub1y subW/figW subH/figH])
plot_NaCl_phase_diagram(Phase)

subplot('position',[sub2x sub2y subW/figW subH/figH])
plotTX_diagram(Phase,-.5,1.5)

subplot('position',[sub3x sub3y subW/figW subH/figH])
HDmax = round(Phase.invHX.HD2l+1);
plot_HX_diagram(Phase,Phase.Xe,-1,HDmax)

subplot('position',[sub4x sub4y subW/figW subH/figH])
CDmax = ceil(max(Phase.invHC.CD3l,Phase.invHC.CD3s)*10)/10;
plot_HC_diagram(Phase,CDmax,-1,HDmax)


