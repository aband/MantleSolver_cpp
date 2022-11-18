function [ru] = weno_reconst(x0,dx,coeff,p,uBarStencil,stencil)

   ru = 0.0;

   r = stencil(2) - stencil(1) + 1;

   for k=1:r
       for i = 1:r
           ru = ru + uBarStencil(k) * coeff(k,i)*((p-x0)/dx).^(i-1);
       end
   end

end

function [sol] = polyn(stencil)

r = stencil(2) - stencil(1) + 1;

M = zeros(r,r);

sol = zeros(r,r);

for k = 1:r
   for j = 1:r
       xleft  = stencil(1)+(j-1)-0.5;
       xright = xleft+1;

       for i=1:r
           M(j,i) = integral(@(x) x.^(i-1), xleft, xright);
       end
   end

   B = zeros(r,1);
   B(k) = 1;
   sol(k,:) = M\B;
end

end

% ==================================

%{
 {function [val3M,val3P] = polyn3(u,iLL,iL,i,iR,iRR)
 {  val3M(1) = (-1/6)*u(iLL) + ( 5/6)*u(iL) + ( 1/3)*u(i);
 {  val3P(1) = ( 1/3)*u(iLL) + (-7/6)*u(iL) + (11/6)*u(i);
 {
 {  val3M(2) = ( 1/3)*u(iL) + (5/6)*u(i) + (-1/6)*u(iR);
 {  val3P(2) = (-1/6)*u(iL) + (5/6)*u(i) + ( 1/3)*u(iR);
 {
 {  val3M(3) = (11/6)*u(i) + (-7/6)*u(iR) + ( 1/3)*u(iRR);
 {  val3P(3) = ( 1/3)*u(i) + ( 5/6)*u(iR) + (-1/6)*u(iRR);
 {end
 {
 {
 {function [val5M,val5P] = polyn5(u,iLL,iL,i,iR,iRR)
 {  val5M = (-1/20)*u(iLL) + (  9/20)*u(iL) + (47/60)*u(i) + (-13/60)*u(iR) + ( 1/30)*u(iRR);
 {  val5P = ( 1/30)*u(iLL) + (-13/60)*u(iL) + (47/60)*u(i) + (  9/20)*u(iR) + (-1/20)*u(iRR);
 {end
 {
 {function [val5M,val5P] = polyn5s(u,iLLL,iLL,iL,i,iR,iRR,iRRR)
 {
 {  val5M(1) = (-1/20)*u(iLL) + (  9/20)*u(iL) + (47/60)*u(i) + (-13/60)*u(iR) + ( 1/30)*u(iRR);
 {  val5P(1) = ( 1/30)*u(iLL) + (-13/60)*u(iL) + (47/60)*u(i) + (  9/20)*u(iR) + (-1/20)*u(iRR);
 {
 {  val5M(2) = (1/5)*u(iL) + (77/60)*u(i) + (-43/60)*u(iR) + (17/60)*u(iRR) + (-1/20)*u(iRRR);
 {  val5P(2) = (-1/20)*u(iL) + (9/20)*u(i) + (47/60)*u(iR) + (-13/60)*u(iRR) + (1/30)*u(iRRR);
 {
 {  val5M(3) = (1/5)*u(iLLL) + (-21/20)*u(iLL) + (137/60)*u(iL) + (-163/60)*u(i) + (137/60)*u(iR);
 {  val5P(3) = (-1/20)*u(iLLL) + (17/60)*u(iLL) + (-43/60)*u(iL) + (77/60)*u(i) + (1/5)*u(iR);
 {end
 {
 {function [val7M,val7P] = polyn7(u,iLLL,iLL,iL,i,iR,iRR,iRRR)
 {  val7M = (1/105)*u(iLLL) + (-19/210)*u(iLL) + (107/210)*u(iL) + (319/420)*u(i) + (-101/420)*u(iR) + (5/84)*u(iRR) + (-1/140)*u(iRRR);
 {  val7P = (-1/140)*u(iLLL) + (5/84)*u(iLL) + (-101/420)*u(iL) + (319/420)*u(i) + (107/210)*u(iR) + (-19/210)*u(iRR) + (1/105)*u(iRRR);
 {end
 {
 {function [sigma3] = smoothness3(u,iLL,iL,i,iR,iRR)
 {  sigma3(1) = 1/2 * ((u(iLL)-u(i))^2/4 + (u(iL)-u(i))^2);
 {  sigma3(2) = 1/2 * ((u(i)-u(iL))^2 + (u(i)-u(iR))^2);
 {  sigma3(3) = 1/2 * ((u(iRR)-u(i))^2/4 + (u(iR)-u(i))^2);
 {end
 {
 {function [sigma5] = smoothness5(u,iLL,iL,i,iR,iRR)
 {  sigma5(1) = 1/4 * (((u(iLL)-u(i))/2)^2 + ((u(iL)-u(i)))^2 + ...
 {                    ((u(iRR)-u(i))/2)^2 + ((u(iR)-u(i)))^2) ;
 {
 {  sigma5(2) = 1/4 * (((u(iRR)-u(i))/2)^2 + ((u(iL)-u(i)))^2 + ...
 {                    ((u(iRRR)-u(i))/3)^2 + ((u(iR)-u(i)))^2) ;
 {
 {  sigma5(3) = 1/4 * (((u(iLLL)-u(i))/3)^2 + ((u(iLL)-u(i))/2)^2 + ...
 {                    ((u(iL)-u(i)))^2 + ((u(iR)-u(i)))^2) ;
 {
 {end
 {
 {function [sigma7] = smoothness7(u,iLLL,iLL,iL,i,iR,iRR,iRRR)
 {  sigma7 = 1/6 * ( ((u(iRRR)-u(i))/3)^2 + ((u(iRR)-u(i))/2)^2 + ((u(iR)-u(i)))^2 + ...
 {                   ((u(iL)-u(i)))^2 + ((u(iLL)-u(i))/2)^2 + ((u(iLLL)-u(i))/3)^2 );
 {
 {end
 %}


