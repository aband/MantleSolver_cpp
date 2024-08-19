function [] = plot_eutectic()

Xe = 0.35; % Eutectic mass fraction

T1 = 2053;
T2 = 1753;

Te = 1560;

Tm = max(T1,T2);

T1 = (T1-Te)/(Tm-Te);
T2 = (T2-Te)/(Tm-Te);
Te = (Te-Te)/(Tm-Te);

% Weight fraction Temperature plot
figure
plot([0 Xe]*100, [T1, Te], 'k-'), hold on
plot([Xe 1]*100, [Te, T2], 'k-')
plot([0  1]*100, [Te, Te], 'k-')
plot([0  Xe]*100, [-0.2, -0.2], 'k-')
plot([0 0 Xe 1 1]*100, [Te T1 Te T2 Te], 'o', 'markerfacecolor', 'w', 'markeredgecolor', 'r')
text(.9, Te+0.04,'1','fontsize',12,'color','r')
text(Xe*100+0.2, Te+0.05,'3','fontsize',12,'color','r')
text(3, T1-0.04,'2','fontsize',12,'color','r')
text(0,-0.1,'I','fontsize',12,'color','b');
text(50,-0.1,'II','fontsize',12,'color','b');
text(3,0.5,'IV','fontsize',12,'color','b');
text(80,0.2,'V','fontsize',12,'color','b');
text(50,0.5,'VI','fontsize',12,'color','b');
text(Xe*100, Te,'III','fontsize',12,'color','b')
hold off

axis square
xlabel 'opx [wt%]'
ylabel 'Dimensionless Temp'

% bulk composition Temperature plot switch to left side of eutectic
figure
plot([0 Xe]/Xe, [T1, Te], 'k-'), hold on
plot([0  Xe]/Xe, [-0.2, -0.2], 'k-')
plot([0 0 Xe]/Xe, [Te T1 Te], 'o', 'markerfacecolor', 'w', 'markeredgecolor', 'r')
plot([0  Xe]/Xe, [Te, Te], 'k-')
text(0.01, Te+0.04,'1','fontsize',12,'color','r')
text(Xe/Xe-0.1, Te+0.05,'3','fontsize',12,'color','r')
text(0.02, T1-0.1,'2','fontsize',12,'color','r')
text(0,-0.1,'I','fontsize',12,'color','b');
text(0.5,-0.1,'II','fontsize',12,'color','b');
text(0.2,0.4,'IV','fontsize',12,'color','b');
text(0.8,0.5,'VI','fontsize',12,'color','b');
text(1, Te,'III','fontsize',12,'color','b')
hold off

axis square
xlabel 'Bulk composition'
ylabel 'Dimensionless Temp'

% H-X plot 
L = 0.6

figure
plot([0,Xe]/Xe,[T1+L,Te+L],'k-'), hold on
plot([0,Xe]/Xe,[Te,Te+L],'k-')
plot([0,Xe]/Xe, [-0.2, -0.2], 'k-')
plot([0,Xe]/Xe,[Te,Te],'k-')
text(0.5,-0.1,'II','fontsize',12,'color','b');
text(0.3,0.6,'IV','fontsize',12,'color','b');
text(0.8,0.3,'III','fontsize',12,'color','b');
text(0.6,1.3,'VI','fontsize',12,'color','b');
hold off

axis square
xlabel 'Bulk composition'
ylabel 'Dimensionless enthalpy'

