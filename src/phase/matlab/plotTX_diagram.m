function [] = plotTX_diagram(Phase,Tmin,Tmax)
col = marc_colors();
% Unpack for readability
Xe = Phase.Xe;

% Phase field boundaries:
plot([0,Xe],[0 0],'k-'), hold on
plot([0 Xe],[1 0],'k-')

% Phase field corners:
plot([0 Xe 0],[0 0 1],'bo','markerfacecolor','w','markeredgecolor',col.red)

xlim([0 Xe])
ylim([Tmin Tmax])

xlabel 'mass fraction: X'
ylabel 'Temperature: T_D'
text(.017,0.075,'1','fontsize',16,'color',col.red)
text(.017,1.075,'2','fontsize',16,'color',col.red)
text(Xe+.017,0.075,'3','fontsize',16,'color',col.red)

pbaspect([.8 1 1])
