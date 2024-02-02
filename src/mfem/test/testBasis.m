clc; clear

% Define basis functions
x_v = [-1,1,1,-1];
y_v = [-1,-1,1,1];

% -- diagnol normal unit vector and linear functions
nu_d1 = [y_v(1)-y_v(3),x_v(3)-x_v(1)]/sqrt((x_v(3)-x_v(1))^2+(y_v(3)-y_v(1))^2);
nu_d2 = [y_v(2)-y_v(4),x_v(4)-x_v(2)]/sqrt((x_v(2)-x_v(4))^2+(y_v(4)-y_v(2))^2);

lam_d1 = @(x) -1*(x-[x_v(1),y_v(1)])*nu_d1';
lam_d2 = @(x) -1*(x-[x_v(2),y_v(2)])*nu_d2';

% -- edge local normal unit vector and linear functions
nu_1 = [y_v(1)-y_v(4),x_v(4)-x_v(1)]/sqrt((y_v(1)-y_v(4))^2+(x_v(1)-x_v(4))^2);
nu_2 = [y_v(2)-y_v(1),x_v(1)-x_v(2)]/sqrt((y_v(1)-y_v(2))^2+(x_v(1)-x_v(2))^2);
nu_3 = [y_v(3)-y_v(2),x_v(2)-x_v(3)]/sqrt((y_v(2)-y_v(3))^2+(x_v(2)-x_v(3))^2);
nu_4 = [y_v(4)-y_v(3),x_v(3)-x_v(4)]/sqrt((y_v(4)-y_v(3))^2+(x_v(3)-x_v(4))^2);

lam_1 = @(x) -(x-[x_v(1),y_v(1)])*nu_1';
lam_2 = @(x) -(x-[x_v(2),y_v(2)])*nu_2';
lam_3 = @(x) -(x-[x_v(3),y_v(3)])*nu_3';
lam_4 = @(x) -(x-[x_v(4),y_v(4)])*nu_4';

vec1 = [x_v(1),y_v(1)];
vec2 = [x_v(2),y_v(2)];
vec3 = [x_v(3),y_v(3)];
vec4 = [x_v(4),y_v(4)];

edge1 = 0.5*(vec1+vec4);
edge2 = 0.5*(vec2+vec1);
edge3 = 0.5*(vec3+vec2);
edge4 = 0.5*(vec4+vec3);

% -- ausillary rational function derivative
r  = @(x) (lam_3(x)*lam_4(x))/(lam_3([x_v(1),y_v(1)])*lam_4([x_v(1),y_v(1)])) - ...
          (lam_1(x)*lam_4(x))/(lam_1([x_v(2),y_v(2)])*lam_4([x_v(2),y_v(2)])) + ...
          (lam_1(x)*lam_2(x))/(lam_1([x_v(3),y_v(3)])*lam_2([x_v(3),y_v(3)])) - ...
          (lam_2(x)*lam_3(x))/(lam_2([x_v(4),y_v(4)])*lam_3([x_v(4),y_v(4)])) ;

dr = @(x) (lam_3(x)*(-nu_4)+lam_4(x)*(-nu_3))/(lam_3(vec1)*lam_4(vec1)) -... 
          (lam_4(x)*(-nu_1)+lam_1(x)*(-nu_4))/(lam_4(vec2)*lam_1(vec2)) +...
          (lam_1(x)*(-nu_2)+lam_2(x)*(-nu_1))/(lam_1(vec3)*lam_2(vec3)) -...
          (lam_2(x)*(-nu_3)+lam_3(x)*(-nu_2))/(lam_2(vec4)*lam_3(vec4));

psi_1 = @(x) lam_2(x).*lam_4(x).*lam_3(x)/(lam_3(x)+lam_1(x));
psi_2 = @(x) lam_3(x).*lam_1(x).*lam_4(x)/(lam_4(x)+lam_2(x));
psi_3 = @(x) lam_4(x).*lam_2(x).*lam_1(x)/(lam_1(x)+lam_3(x));
psi_4 = @(x) lam_1(x).*lam_3(x).*lam_2(x)/(lam_2(x)+lam_4(x));

phi_e11 = @(x) psi_1(x)/psi_1(edge1);
phi_e21 = @(x) psi_2(x)/psi_2(edge2);
phi_e31 = @(x) psi_3(x)/psi_3(edge3);
phi_e41 = @(x) psi_4(x)/psi_4(edge4);

R = @(x) r(x) - r(edge1)*phi_e11(x)-...
                r(edge2)*phi_e21(x)-...
                r(edge3)*phi_e31(x)-...
                r(edge4)*phi_e41(x);

phi_v1 = @(x) (lam_d2(x)-0.5*lam_d2(vec3)*(1+R(x)))/(lam_d2(vec1)-lam_d2(vec3));
phi_v2 = @(x) (lam_d1(x)-0.5*lam_d1(vec4)*(1-R(x)))/(lam_d1(vec2)-lam_d1(vec4));
phi_v3 = @(x) (lam_d2(x)-0.5*lam_d2(vec1)*(1+R(x)))/(lam_d2(vec3)-lam_d2(vec1));
phi_v4 = @(x) (lam_d1(x)-0.5*lam_d1(vec2)*(1-R(x)))/(lam_d1(vec4)-lam_d1(vec2));

% ----

dphi_11x = @(x) (-(lam_1(x)+lam_3(x))*(lam_2(x)*lam_4(x)*nu_3(1) + lam_2(x)*lam_3(x)*nu_4(1) + lam_3(x)*lam_4(x)*nu_2(1)) + lam_2(x)*lam_3(x)*lam_4(x)*(nu_1(1)+nu_3(1))) /...
               ((lam_1(x)+lam_3(x))^2*psi_1(edge1));

dphi_21x = @(x) (-(lam_2(x)+lam_4(x))*(lam_1(x)*lam_4(x)*nu_3(1) + lam_1(x)*lam_3(x)*nu_4(1) + lam_3(x)*lam_4(x)*nu_1(1)) + lam_4(x)*lam_3(x)*lam_1(x)*(nu_2(1)+nu_4(1))) /...
               ((lam_2(x)+lam_4(x))^2*psi_2(edge2));

dphi_31x = @(x) (-(lam_1(x)+lam_3(x))*(lam_2(x)*lam_4(x)*nu_1(1) + lam_2(x)*lam_1(x)*nu_4(1) + lam_1(x)*lam_4(x)*nu_2(1)) + lam_2(x)*lam_1(x)*lam_4(x)*(nu_1(1)+nu_3(1))) /...
               ((lam_1(x)+lam_3(x))^2*psi_3(edge3));

dphi_41x = @(x) (-(lam_2(x)+lam_4(x))*(lam_2(x)*lam_1(x)*nu_3(1) + lam_2(x)*lam_3(x)*nu_1(1) + lam_3(x)*lam_1(x)*nu_2(1)) + lam_2(x)*lam_3(x)*lam_1(x)*(nu_2(1)+nu_4(1))) /...
               ((lam_2(x)+lam_4(x))^2*psi_4(edge4)); 

% ---------------------

dphi_11y = @(y) (-(lam_1(y)+lam_3(y))*(lam_2(y)*lam_4(y)*nu_3(2) + lam_2(y)*lam_3(y)*nu_4(2) + lam_3(y)*lam_4(y)*nu_2(2)) + lam_2(y)*lam_3(y)*lam_4(y)*(nu_1(2)+nu_3(2))) /...
               ((lam_1(y)+lam_3(y))^2*psi_1(edge1));

dphi_21y = @(y) (-(lam_2(y)+lam_4(y))*(lam_1(y)*lam_4(y)*nu_3(2) + lam_1(y)*lam_3(y)*nu_4(2) + lam_3(y)*lam_4(y)*nu_1(2)) + lam_4(y)*lam_3(y)*lam_1(y)*(nu_2(2)+nu_4(2))) /...
               ((lam_2(y)+lam_4(y))^2*psi_2(edge2));

dphi_31y = @(y) (-(lam_1(y)+lam_3(y))*(lam_2(y)*lam_4(y)*nu_1(2) + lam_2(y)*lam_1(y)*nu_4(2) + lam_1(y)*lam_4(y)*nu_2(2)) + lam_2(y)*lam_1(y)*lam_4(y)*(nu_1(2)+nu_3(2))) /...
               ((lam_1(y)+lam_3(y))^2*psi_3(edge3));

dphi_41y = @(y) (-(lam_2(y)+lam_4(y))*(lam_2(y)*lam_1(y)*nu_3(2) + lam_2(y)*lam_3(y)*nu_1(2) + lam_3(y)*lam_1(y)*nu_2(2)) + lam_2(y)*lam_3(y)*lam_1(y)*(nu_2(2)+nu_4(2))) /...
               ((lam_2(y)+lam_4(y))^2*psi_4(edge4)); 


dR = @(x) dr(x) - r(edge1) * [dphi_11x(x),dphi_11y(x)] - ...
                  r(edge2) * [dphi_21x(x),dphi_21y(x)] - ...
                  r(edge3) * [dphi_31x(x),dphi_31y(x)] - ...
                  r(edge4) * [dphi_41x(x),dphi_41y(x)] ;

% -- derivative of basis

grad1 = @(vec) 1.0/(lam_d2([x_v(1),y_v(1)])-lam_d2([x_v(3),y_v(3)])) * ((-1*nu_d2) - 0.5*lam_d2([x_v(3),y_v(3)])* dR(vec));
grad2 = @(vec) 1.0/(lam_d1([x_v(2),y_v(2)])-lam_d1([x_v(4),y_v(4)])) * ((-1*nu_d1) + 0.5*lam_d1([x_v(4),y_v(4)])* dR(vec));
grad3 = @(vec) 1.0/(lam_d2([x_v(3),y_v(3)])-lam_d2([x_v(1),y_v(1)])) * ((-1*nu_d2) - 0.5*lam_d2([x_v(1),y_v(1)])* dR(vec));
grad4 = @(vec) 1.0/(lam_d1([x_v(4),y_v(4)])-lam_d1([x_v(2),y_v(2)])) * ((-1*nu_d1) + 0.5*lam_d1([x_v(2),y_v(2)])* dR(vec));

grad_bubble_1x = @(x) -[dphi_11x(x)*nu_1(1),dphi_11y(x)*nu_1(1)];
grad_bubble_1y = @(x) -[dphi_11x(x)*nu_1(2),dphi_11y(x)*nu_1(2)];
grad_bubble_2x = @(x) -[dphi_21x(x)*nu_2(1),dphi_21y(x)*nu_2(1)];
grad_bubble_2y = @(x) -[dphi_21x(x)*nu_2(2),dphi_21y(x)*nu_2(2)];
grad_bubble_3x = @(x) [dphi_31x(x)*nu_3(1),dphi_31y(x)*nu_3(1)];
grad_bubble_3y = @(x) [dphi_31x(x)*nu_3(2),dphi_31y(x)*nu_3(2)];
grad_bubble_4x = @(x) [dphi_41x(x)*nu_4(1),dphi_41y(x)*nu_4(1)];
grad_bubble_4y = @(x) [dphi_41x(x)*nu_4(2),dphi_41y(x)*nu_4(2)];

% Test of all the functions
seed   = 5;
startp = -1;
endp   = 1;

edge   = 0;
shift  = 0;

dof    = edge + 4*shift;

dh = (endp - startp)/seed;

for jj = seed : -1 : 0
yy = startp + jj*dh;
for ii = 0:1:seed
xx = startp + ii*dh; 
currentVec = [xx,yy];

gradVec1 = grad_bubble_4x(currentVec);
gradVec2 = grad_bubble_4y(currentVec);

fprintf('(%.2f , %.2f, %.2f, %.2f) ', gradVec1(1), gradVec1(2),...
                                      gradVec2(1), gradVec2(2));
end
fprintf('\n');
end
