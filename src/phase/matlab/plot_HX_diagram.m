function []=plot_HX_diagram(Phase,Xmax,HDmin,HDmax)
col = marc_colors();

Xe   = Phase.Xe;
HD2l = Phase.invHX.HD2l;
HD3l = Phase.invHX.HD3l;
HDe  = Phase.invHX.HDe;
HDl  = Phase.invHX.HDl;

X = linspace(0,Xe,10);

plot([0 Xe],[0 0],'k-'), hold on
plot(Xe*[1 1],[0 HD3l],'k-')
plot(X,HDe(X),'k')
plot(X,HDl(X),'k')

% Corners
plot(Xe*[1 1 0 0 0],[0 HD3l 0 HD2l 1],'o','markerfacecolor','w','markeredgecolor',col.red)

text(.017,.4,'1','fontsize',16,'color',col.red)
text(.017,HD2l+0.4,'2l','fontsize',16,'color',col.red)
text(.017,1.4,'2s','fontsize',16,'color',col.red)
text(Xe+.017,0.4,'3s','fontsize',16,'color',col.red)
text(Xe+.017,HD3l+0.4,'3l','fontsize',16,'color',col.red)


xlim([0 Xmax])
ylim([HDmin HDmax])
xlabel 'mass fraction: X'
ylabel 'enthalpy: H_D'
pbaspect([.8 1 1])