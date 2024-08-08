function [] = plot_eutectic()

Xe = 0.3; % Eutectic mass fraction

T1 = 1333;
T2 = 1350;

Te = 1280;

Tm = max(T1,T2);

T1 = (T1-Te)/(Tm-Te);
T2 = (T2-Te)/(Tm-Te);
Te = (Te-Te)/(Tm-Te);

plot([0 Xe]*100, [T1, Te], 'k-'), hold on
plot([Xe 1]*100, [Te, T2], 'k-')
plot([0  1]*100, [Te, Te], 'k-')
plot([0  0]*100, [-0.03, -0.03], 'k-')
plot([0 0 Xe 1 1]*100, [Te T1 Te T2 Te], 'o', 'markerfacecolor', 'w', 'markeredgecolor', 'r')
text(.9, Te+0.04,'1','fontsize',12,'color','r')

axis square
xlabel '[wt%]'
ylabel 'Dimensionless Temp'
