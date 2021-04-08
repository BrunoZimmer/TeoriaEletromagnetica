% Um disco de raio a pertence ao plano xy, com o eixo z passando pelo seu centro. 
% Uma carga superficial de densidade uniforme ?s está presente no disco, que gira 
% em volta do eixo z em uma velocidade angular de w rad/s. Calcule numericamente 
% o potencial vetorial magnético A e apartir dele a densidade de fluxo mangnético B, 
% em qualquer posição.  Faça uma representação gráfica de A e de B no plano xz.
% 
% Afim de valiar sua resposta calcule a magnitude da intensidade de campo 
% magnético  na posição x = 0, y = 0, z= 6,9 m. Considere a = 6,1 m, 
% ?s = 9,4 C/m2 e ? = 7,1 rad/s

clc
clear all
close all

%% Variáveis iniciais

u0 = 4*pi*10^-7; % permeablidade do espaço livre (em H/m)
raio = 6.1; % raio do anel
w = 7.1; %velocidade angular [rad/s]
ps = 9.4;%densidade superficial [C/m^2]
zmax = 6.9;

%% Espaços de calculos

D = 1.5*raio; % limite para as coordenadas x, y e z
passo = zmax/25; % passo

dx = passo;
dy = passo;
dz = passo;

x = -D:dx:D; %variacao da coordenada x 
y = -D:dy:D; %variacao da coordenada y 
z = -2*zmax:dz:2*zmax; %variacao da coordenada z

dphi = 2*pi/14;
drho = raio/20;
phil = 0:dphi:2*pi-dphi/2; % notar que excluímos 2*pi, pois já¡ temos 0.
rhol = 0:drho:raio;

%% Descobrindo valor de A

%zerando o potencial vetorial magnético
Ax = zeros(length(x),length(y),length(z));
Ay = zeros(length(x),length(y),length(z));
Az = zeros(length(x),length(y),length(z));

for i = 1:length(x)% varre a coordenada x onde H será¡ calculado
    disp(i)
    for j = 1: length(y)  % varre a coordenada y onde H será¡ calculado
        for k = 1:length(z) % varre a coordenada z onde H será¡ calculado
            
            for m = 1:length(phil)  % varre a coordenada phil, onde existe a corrente
                for n =1:length(rhol)
                    % vetor posição onde calculamos A
                    r = [x(i),y(j),z(k)];
                    
                    % rl aponta para o diferencial dL 
                    [rl(1), rl(2), rl(3)] = pol2cart(phil(m), rhol(n), 0);
                    
                    %[K(1), K(2), K(3)] = pol2cart(pi/2+phil(m), rhol(n)*ps*w, 0); %â_phi = â_y*cos(phi) - â_x*sen(phi)
                    %pol2cart confuso
                    K = [rhol(n)*ps*w*(-sin(phil(m))), rhol(n)*ps*w*cos(phil(m)), 0];
                    
                    dS = rhol(n)*drho*dphi;
                    KdS = K*dS;
                    % pequena parte para campo A de K em cada posição
                    dA = KdS/sqrt( (r-rl)*(r-rl)');
                    
                    Ax(i,j,k)=Ax(i,j,k)+dA(1);%soma em Ax a componente x da contribuição em dA
                    Ay(i,j,k)=Ay(i,j,k)+dA(2);%soma em Ay a componente y da contribuição em dA
                    Az(i,j,k)=Az(i,j,k)+dA(3);%soma em Az a componente z da contribuição em dA
                end
            end
        end
    end
end

Ax=(u0/(4*pi))*Ax;
Ay=(u0/(4*pi))*Ay;
Az=(u0/(4*pi))*Az;

%% Graficos de A

xzero=ceil(length(x)/2);
yzero=ceil(length(y)/2);
zzero=ceil(length(z)/2);

%Gráfico XZ de A em y = 0
figure(1);

[X,Z] = meshgrid(x,z);
quiver(X, Z, squeeze(Ax(:,yzero,:))', squeeze(Az(:,yzero,:))');  
xlabel('eixo x (m)')
ylabel('eixo z (m)')
title('A, plano XZ (y=0)');

%Gráfico YZ de A em x = 0
figure(2);

[Y,Z] = meshgrid(y,z);
quiver(Y, Z, squeeze(Ay(xzero,:,:))', squeeze(Az(xzero,:,:))');
xlabel('eixo y (m)')
ylabel('eixo z (m)')
title('A, plano YZ (x=0)');

%Gráfico XY de A em z = 0
figure(3);

[X,Y] = meshgrid(x,y);
quiver(X, Y, squeeze(Ax(:,:,zzero))', squeeze(Ay(:,:,zzero))');  %faz o gráfico de A no plano y = 0
xlabel('eixo x (m)')
ylabel('eixo y (m)')
title('A, plano XY (z=0)');


%% Descobrindo valor de B

%Derivada é um espaço menor tirando o "cantos" do espaço
xn = (x(2:end-1));
yn = (y(2:end-1));
zn = (z(2:end-1));

%Zerando Densisadede de fluxo magnetico(B)
Bx = zeros (length(xn),length(yn),length(zn));
By = zeros (length(xn),length(yn),length(zn));
Bz = zeros (length(xn),length(yn),length(zn));

%B a partir do rotacional de A 
%rot A = (dAz/dy-dAy/dz) ax + (dAx/dz-dAz/dx) ay + (dAy/dx-dAx/dy) az

for i = 2:length(x)-1% varre a coordenada x onde B será calculado
    disp(i)
    for j = 2 :length(y)-1 % varre a coordenada y onde B será calculado
        for k = 2: length(z)-1  % varre a coordenada z onde E será calculado
            
            dAz_dy= (Az(i,j+1,k)-Az(i,j-1,k))/2/dy; %derivada de Az em relação a y
            dAy_dz = (Ay(i,j,k+1)-Ay(i,j,k-1))/2/dz; %derivada de Ay em relação a z
            
            dAx_dz=(Ax(i,j,k+1)-Ax(i,j,k-1))/2/dz; %derivada de Ax em relação a z
            dAz_dx=(Az(i+1,j,k)-Az(i-1,j,k))/2/dx;%derivada de Az em relação a x
            
            dAy_dx=(Ay(i+1,j,k)-Ay(i-1,j,k))/2/dx; %derivada de Ay em relação a x
            dAx_dy=(Ax(i,j+1,k)-Ax(i,j-1,k))/2/dy;%derivada de Ax em relação a y
            
            Bx(i-1,j-1,k-1) = dAz_dy-dAy_dz; %a componente x de B é dAz/dy-dAy/dz
            By(i-1,j-1,k-1) = dAx_dz-dAz_dx;%a componente y de B é dAx/dz-dAz/dx
            Bz(i-1,j-1,k-1) = dAy_dx-dAx_dy;%a componente z de B é dAy/dx-dAx/dy
        end
    end
end

%% Graficos B

xzero2 = ceil(length(xn)/2);
yzero2 = ceil(length(yn)/2);
zzero2 = ceil(length(zn)/2);

%Gráfico XZ de B em y = 0
figure(4)
[X,Z] = meshgrid(xn,zn);
quiver(X,Z,squeeze(Bx(:,yzero2,:))' , squeeze(Bz(:,yzero2,:))')
xlabel('eixo x (m)')
ylabel('eixo Z (m)')
title('B, plano XZ (y=0)');

%% Resultado no ponto
%A no ponto (0; 0; 6.9)

Hresp = 0;
Hresp = Hresp + Bx(xzero2,yzero2,zzero2*1.5+1);
Hresp = Hresp + By(xzero2,yzero2,zzero2*1.5+1);
Hresp = Hresp + Bz(xzero2,yzero2,zzero2*1.5+1);

Hresp = Hresp/u0;

disp(Hresp);%-0,93086 A/m.

