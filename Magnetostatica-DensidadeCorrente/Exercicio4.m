%Seja J = 25/p ap -20/(p2 + 0,01) az A/m2. Determine numericamente o campo
%vetorial J para uma região do espaço que englobe o plano z = 9,4 e
%0 < p < 2*3,0. Trace o gráfico de J no plano yz (phi=pi/2 ou x=0). Trace o gráfico 
%de J no plano xy (z=0). Trace o gráfico de J no plano z = 9,4. Trace o 
%gráfico de J no plano z = -9,4.Com base nos dados calculados, integre 
%numericamente o fluxo de J através do plano z = 9,4, na região em que
%p < 3,0 m, determinando a corrente que passa por essa região.

clc
clear all
close all

%% Condicoes iniciais
u0=4*pi*10^(-7); % Permeabilidade magnetica do espaÃ§o livre (em H/m)

zm = 9.4; %magnitude da analise em z
pm = 2*(3);%magnitude da analise em p

passo=pm/120; % passo

dx=passo;
dy=passo;
dz=passo;

x= -pm:dx:pm; %variacao da coordenada x
y= -pm:dy:pm; %variacao da coordenada y
z= -zm:dz:zm; %variacao da coordenada z


% Inicializa a densidade volumetrica de corrente
Jx(:,:,:) = zeros(length(x),length(y),length(z)); 
Jy(:,:,:) = zeros(length(x),length(y),length(z));
Jz(:,:,:) = zeros(length(x),length(y),length(z));


% Ã?ndice do valor mÃ©dio dos vetores de x, y e z.
% (Neste exemplo, sÃ£o os Ã­ndices em que x=y=z=0.)

xmedio=ceil(length(x)/2);
ymedio=ceil(length(y)/2);
zmedio=ceil(length(z)/2);

%% Calcular J em todo o espaço xyz
for i = 1:length(x)% varre a coordenada x onde B serÃ¡ calculado
%   disp(i)
  for j = 1: length(y)  % varre a coordenada y onde B serÃ¡ calculado
    for k = 1:length(z) % varre a coordenada z onde B serÃ¡ calculado
   
      %J = 25/p ap -20/(p2 + 0,01) az A/m2
      %J = 25/(x^2+y^2)^(1/2) ap -20/((x^2+y^2) + 0,01) az A/m2
      
      %ax = ap*cos(atan(y/x))
      % cos(atan2( y(j), x(i))) ou x(i)/sqrt(x(i)^2 + y(j)^2)
      % sin(atan2( y(j), x(i))) ou y(j)/sqrt(x(i)^2 + y(j)^2)
      Jx(i,j,k) = (25/(x(i)^2+y(j)^2)^(1/2)) * x(i)/sqrt(x(i)^2 + y(j)^2); %x(i)/sqrt(x(i)^2 + y(j)^2)
      Jy(i,j,k) = (25/(x(i)^2+y(j)^2)^(1/2)) * y(j)/sqrt(x(i)^2 + y(j)^2);
      Jz(i,j,k) = (-20)/((x(i)^2+y(j)^2) + 0.01);

    end
  end
end

%% Graficos J
%gráfico de J no plano yz (phi=pi/2)

% phi = atan(y/x)
% y = x * tan(phi)

figure(1);

[Y,Z] = meshgrid(y,z);
quiver(Y, Z, squeeze(Jy(xmedio,:,:))', squeeze(Jz(xmedio,:,:))'); 
hold on
contour(y, z, squeeze(Jx(xmedio,:,:))', 20);

xlabel('eixo y (m)')
ylabel('eixo z (m)')
title('J no plano yz (phi=pi/2)');

%gráfico de J no plano xy (z=0)

figure(2);

[X,Y] = meshgrid(x,y);
quiver(X, Y, squeeze(Jx(:,:,zmedio))', squeeze(Jy(:,:,zmedio))');
hold on
contour(x, y, squeeze(Jz(:,:,zmedio))', 20);

xlabel('eixo x (m)')
ylabel('eixo y (m)')
title('J no plano xy (z=0)');

%gráfico de J no plano z = 9,4

figure(3);

[X,Y] = meshgrid(x,y);
quiver(X, Y, squeeze(Jx(:,:,length(z)))', squeeze(Jy(:,:,length(z)))'); 
hold on
contour(x, y, squeeze(Jz(:,:,length(z)))', 20);

xlabel('eixo x (m)')
ylabel('eixo y (m)')
title('J no plano z = 9,4');

%gráfico de J no plano z = -9,4

figure(4);

[X,Y] = meshgrid(x,y);
quiver(X, Y, squeeze(Jx(:,:,1))', squeeze(Jy(:,:,1))'); 
hold on
contour(x, y, squeeze(Jz(:,:,1))', 20);

xlabel('eixo x (m)')
ylabel('eixo y (m)')
title('J no plano z = -9,4');

%% Calculando corrente 
%corrente em z = 9,4 e 0 < p < 3,0 m

% Inicializa a corrente
cur = 0;

for i = 2:length(x)% varre a coordenada x
%   disp(i)
  for j = 2: length(y)% varre a coordenada y
      if ( ((x(i)^2+y(j)^2)^(1/2)) <= 3 && ((x(i)^2+y(j)^2)^(1/2)) > 0  )
        cur = cur+Jx(i,j,length(z))*dx*dy;
        cur = cur+Jy(i,j,length(z))*dx*dy;
        cur = cur+Jz(i,j,length(z))*dx*dy;
      end  
      
  end
end
disp(cur);%-427,48
