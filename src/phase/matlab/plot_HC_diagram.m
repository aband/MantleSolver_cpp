function [] = plot_HC_diagram(Phase,CDmax,HDmin,HDmax)
% col = marc_colors();

HDe  = Phase.invHC.HDe;
HDl  = Phase.invHC.HDl;
HDb  = Phase.invHC.HDb;
CD3s = Phase.invHC.CD3s;
CD3l = Phase.invHC.CD3l;
HD2l = Phase.invHC.HD2l;
HD3l = Phase.invHC.HD3l;

CD  = linspace(0,CD3l,10);
CDb = linspace(CD3s,CD3l,10);

% Boundaries
plot([0 CD3s],[0 0],'k-'), hold on
plot(CD,HDe(CD),'k-')
plot(CD,HDl(CD),'k-')
plot(CDb,HDb(CDb),'k-')

% Corners
plot([0 CD3s CD3l 0 0 ],[0 0 HD3l HD2l 1],'o','markerfacecolor','w','markeredgecolor','r')

text(.03,.4,'1','fontsize',16,'color','r')
text(.03,HD2l+0.4,'2l','fontsize',16,'color','r')
text(.03,1.4,'2s','fontsize',16,'color','r')
text(CD3s+.03,0.4,'3s','fontsize',16,'color','r')
text(CD3l+.03,HD3l+0.4,'3l','fontsize',16,'color','r')
xlim([0 CDmax])
ylim([HDmin HDmax])

xlabel 'composition: C_D'
ylabel 'enthalpy: H_D'
pbaspect([.8 1 1])
