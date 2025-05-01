function [] = plot_eutectic()
col = marc_colors();

Xe = 0.25; % Eutectic mass fraction

T1 = 2000;
T2 = 1800;

Te = 1480;

Tm = max(T1,T2);

T1 = (T1-Te)/(Tm-Te);
T2 = (T2-Te)/(Tm-Te);
Te = (Te-Te)/(Tm-Te);

fontSize = 14;

% Weight fraction Temperature plot
figure
plot([0 Xe]*100, [T1, Te], 'k-'), hold on
plot([Xe 1]*100, [Te, T2], 'k-')
plot([0  1]*100, [Te, Te], 'k-')
plot([0  Xe]*100, [-0.2, -0.2], 'k-')
plot([0 0 Xe 1 1, Xe]*100, [Te T1 Te T2 Te,  -0.2], 'o', 'markerfacecolor', 'w', 'markeredgecolor', 'r')
text(.9, Te+0.04,'1','fontsize',12,'color','r')
text(0.03*Xe, Te-0.05,'$T_e$','fontsize',12,'color','r', 'interpreter','latex')
text(Xe*100+0.2, Te+0.05,'3','fontsize',12,'color','r')

text(Xe*100, -0.2+0.05,'$C_e$','fontsize',12,'color','r','interpreter','latex')

text(3, T1-0.04,'2','fontsize',12,'color','r')

text(0,-0.1,'I','fontsize',12,'color','b');
text(50,-0.1,'II','fontsize',12,'color','b');
text(3,0.5,'IV','fontsize',12,'color','b');
text(80,0.2,'V','fontsize',12,'color','b');
text(50,0.5,'VI','fontsize',12,'color','b');
text(Xe*100, Te,'III','fontsize',12,'color','b')
hold off

axis square
xlabel ('opx [wt%]', 'FontSize', fontSize)
ylabel ('Dimensionless Temperature', 'FontSize', fontSize)

%{
% bulk composition Temperature plot switch to left side of eutectic
figure
plot([0 Xe], [T1, Te], 'k-'), hold on
plot([0  Xe], [-0.2, -0.2], 'k-')
plot([0 0 Xe], [Te T1 Te], 'o', 'markerfacecolor', 'w', 'markeredgecolor', 'r')
plot([0  Xe], [Te, Te], 'k-')
text(0.01*Xe, Te+0.04,'1','fontsize',12,'color','r')
text(Xe-0.1*Xe, Te+0.05,'3','fontsize',12,'color','r')
text(0.02*Xe, T1-0.1,'2','fontsize',12,'color','r')
text(0,-0.1,'I','fontsize',12,'color','b');
text(0.5*Xe,-0.1,'II','fontsize',12,'color','b');
text(0.2*Xe,0.4,'IV','fontsize',12,'color','b');
text(0.8*Xe,0.5,'VI','fontsize',12,'color','b');
text(Xe, Te,'III','fontsize',12,'color','b')

% illustration of melting path
plot([0.17 0.17], [-0.2, 0], 'color', col.green, 'linewidth', 2)
text(0.175, -0.16, 'a', 'fontsize', 12, 'color', 'r')
plot([0.17 0.17], [-0.2, 0], 'color', col.green, 'linewidth', 2)
text(0.175, 0.03, 'b', 'fontsize', 12, 'color', 'r')
plot([0.17, Xe], [0,0], 'color', col.green, 'linewidth', 2)
plot([Xe, 0.17], [0,0.32], 'color', col.green, 'linewidth', 2)
plot([0.17, 0.17], [0.32, 1], 'color', col.green, 'linewidth', 2)
text(0.171, 0.28, 'c', 'fontsize', 12, 'color', 'r')
hold off

axis square
xlabel ('Bulk composition', 'FontSize', fontSize)
ylabel ('Dimensionless Temperature', 'FontSize', fontSize)
%}

% H-X plot 
%L = 0.6;
L = 5*10^5 / 1200 / (2000-1480)

%{
figure
plot([0,Xe],[T1+L,Te+L],'k-'), hold on
plot([0,Xe],[Te,Te+L],'k-')
plot([0,Xe], [-0.5, -0.5], 'k-')
plot([0,Xe],[Te,Te],'k-')
text(0.5*Xe,-0.25,'II','fontsize',12,'color','b');
text(0.3*Xe,L,'IV','fontsize',12,'color','b');
text(0.8*Xe,0.3,'III','fontsize',12,'color','b');
text(0.01*Xe,-0.25,'I','fontsize',12,'color','b');
text(0.01*Xe, Te+0.05,'1','fontsize',12,'color','r')
text(0.01*Xe, T1+L+0.05,'2','fontsize',12,'color','r')
text(Xe - 0.03*Xe, Te+0.05,'3','fontsize',12,'color','r')
text(Xe - 0.03*Xe, Te+L+0.05,'4','fontsize',12,'color','r')
plot([0 0 Xe Xe], [Te T1+L Te Te+L], 'o', 'markerfacecolor', 'w', 'markeredgecolor', 'r')
%text(L,1.3,'VI','fontsize',12,'color','b');
%text(0.5*Xe,0.0,'a','fontsize',12,'color','r');%'#A2142F');
%text(0.5*Xe,0.3,'b','fontsize',12,'color','r');%'#A2142F');
%text(0.5*Xe,1.2,'c','fontsize',12,'color','r');%'#A2142F');
hold off

axis square
xlabel ('Bulk composition', 'FontSize', fontSize)
ylabel ('Bulk enthalpy', 'FontSize', fontSize)
%}

%{
% H-C plot
myfunc = @(x) x.*(1-x) + L*x;
xx = linspace(0,1,100);

figure
plot([0,Xe]/Xe,[T1+L,Te+L],'k-'), hold on
plot([0,Xe]/Xe, [-0.2, -0.2], 'k-')
plot([0,Xe]/Xe,[Te,Te],'k-')
%plot(xx,myfunc(xx),'k-');

text(0.5,-0.1,'II','fontsize',12,'color','b');
text(0.3,L,'IV','fontsize',12,'color','b');
text(0.8,0.3,'III','fontsize',12,'color','b');
text(L,1.3,'VI','fontsize',12,'color','b');
%text(0.5,0.55,'b','fontsize',12,'color','#A2142F');
plot([0,Xe]/Xe,[Te,Te+L],'k-')
text(0.5,0.0,'a','fontsize',12,'color','#A2142F');
text(0.5,1.1,'c','fontsize',12,'color','#A2142F');
text(0.5,0.3,'b','fontsize',12,'color','#A2142F');
hold off

axis square
xlabel 'Concentration'
ylabel 'Dimensionless enthalpy'
%}
% 3D plot including pressure
%figure
%[X,p] = meshgrid(0:0.1:1,0:0.2:2);
%
%lp = @(p) 0.6*p
%
%surf1 = @(X,p) lp(p) 
%Z1 = surf1(X,p)
%
%surf2 = @(X,p) (T1+(Te-T1)*X + L)+lp(p)
%Z2 = surf2(X,p)
%
%surf3 = @(X,p) (Te + L*X) + lp(p)
%Z3 = surf3(X,p)
%
%mesh(X,p,Z1), hold on
%mesh(X,p,Z2)
%mesh(X,p,Z3)
%
%axis square
%xlabel 'Concentration'
%ylabel 'Pressure'
%zlabel 'Enthalpy'
