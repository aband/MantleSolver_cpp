function [] = plot_NaCl_phase_diagram(Phase)
col = marc_colors();
Xe = 0.233;
Xs = 0.6;
T1 = Phase.T1;
Te = Phase.Te;
plot([0 Xe]*100,[T1,Te],'k-'), hold on
plot([0 Xs]*100,[Te,Te],'k-')
plot([Xe Xs]*100,[Te T1+75],'k-')
%plot([Xs .3]*100,[T1+0.1 T1+0.1],'k-')
%plot([Xs 0.264]*100,[T1+0.1 1275],'k-')
plot([0 0 Xe]*100,[Te T1 Te],'o','markerfacecolor','w','markeredgecolor',col.red)

ylim([1100 1475])
text(0.15*100,265,'Brine','fontsize',16)
text(0.05*100,257,'Ice + Brine','fontsize',16)
text(0.07*100,248,'Ice + NaCl * 2 H_2O','fontsize',16)
text(.9,Te+1.15,'1','fontsize',16,'color',col.red)
text(.9,T1+1.15,'2','fontsize',16,'color',col.red)
text(100*Xe+.9,Te+1.15,'3','fontsize',16,'color',col.red)


axis square
xlabel 'Olivine [wt%]'
ylabel 'Temperature [K]'
pbaspect([.8 1 1])
