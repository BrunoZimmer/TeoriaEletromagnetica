%TESTES 

%GRAFICO 3D COM FUNÇÃO PRONTA

% x = linspace(-5,5);
% [X,Y] = meshgrid(x);
% M = -log(X.^2+Y.^2);
% %M = exp(-(X.^2+Y.^2)/10);
% figure(1)
% surf(M)
% grid on
% colormap(jet(20))
% colorbar

%GRAFICO 3D COM FUNÇÃO EM DEGRAU

% xd = linspace(-0.01,0.01);
% yd = linspace(-0.01,0.01);
% [x,y] = meshgrid(xd,yd);
% 
% %variaveis do problema
% e0 = 8.854*10^-12;
% pv = 8.1*10^-12;
% b = 7.3*10^-3;
% a = 2.7*10^-3;
% 
% r = sqrt(x.^2+y.^2);
% 
% V = (-(pv*r)/(2*e0))*(0.5 + log(b/a));
% 
% 
% 
% E1 = ((pv.*r)/(2*e0));
% E2 =  ((pv*a^2)/(2*e0)).*(1./r);
% 
% f2 = E1 + (E2-E1).*((a <= r) & (r < b)) -( E1).*(b <= r);
% 
% figure(1)
% surf(x,y,f2)
% grid on
% colormap(jet(20))
% colorbar

% V1 = -((pv.*r^2)/(4*e0));
% V2 =  -((pv*a^2)/(2*e0)).*log(r);

% f1 = V1 + (V2-V1).*((a <= r) & (r < b)) -( V1).*(b <= r);
% 
% figure(1)
% surf(x,y,f1)
% grid on
% colormap(jet(20))
% colorbar

SPHERE3D(1,0,pi,0,pi,1)

