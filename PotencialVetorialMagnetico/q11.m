% Uma lâmina de corrente K flui na região  -a < y < a no plano z =0. Calcule 
% numericamente o potencial vetorial magnético A e a partir desse a densidade 
% de fluxo magnético B, em qualquer posição do espaço. Fação uma representação 
% gráfica de A e de B no plano xz
% Afim de avaliar sua resposta calcule a componente y da intensidade de campo 
% magnético na posição (x = 0, y = 0, z = 7,3 m), considere K = 7,5 ax A/m e 
% a = 3 m.

clc
clear all
close all

%% Variáveis iniciais

u0 = 4*pi*10^-7; % permeablidade do espaço livre (em H/m)
a = 3; % limite da lâmina
k1 = 7.5; %na direção ax
zmax = 7.3;

%% Espaços de calculos

D=2*a; % limite para as coordenadas x, y e z
passo=zmax/20; % passo

dx=passo;
dy=passo;
dz=passo;

x= -D:dx:D; %variacao da coordenada x 
y= -D:dy:D; %variacao da coordenada y
z= -2*zmax:dz:2*zmax; %variacao da coordenada z

xl= -2*D:dx:2*D; %variacao da coordenada x onde está a carga
yl= -a:dy:a; %variacao da coordenada y onde está a carga

%% Descobrindo valor de A

%zerando o potencial vetorial magnético
Ax = zeros(length(x),length(y),length(z));
Ay = zeros(length(x),length(y),length(z));
Az = zeros(length(x),length(y),length(z));

for i = 1:length(x)% varre a coordenada x onde H será¡ calculado
    disp(i)
    for j = 1: length(y)  % varre a coordenada y onde H será calculado
        for k = 1:length(z) % varre a coordenada z onde H será calculado
            
            for m = 1:length(xl)  % varre a coordenada xl, onde existe a corrente
                for n = 1:length(yl)
                    
                    %vetor posição do ponto no espaço
                    r = [x(i),y(j),z(k)];  
                    %vetor posição do ponto na carga
                    rl = [xl(m), yl(n), 0];
                    
                    K = k1;%ax
                    dS = dx*dy;
                    K_dS = K*dS*[1,0,0];
                    
                    % pequena contribuicao para campo A de K em cada posição
                    dA = K_dS/sqrt( (r-rl)*(r-rl)');
                    
                    Ax(i,j,k)=Ax(i,j,k)+dA(1); %soma em Ax a componente x da contribuição em dA
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

xzero=ceil(length(x)/2);
yzero=ceil(length(y)/2);
zzero=ceil(length(z)/2);

%% Graficos A

%Gráfico XZ de A em y = 0
figure(1);
[X,Z] = meshgrid(x,z);
quiver(X, Z, squeeze(Ax(:,yzero,:))', squeeze(Az(:,yzero,:))');
xlabel('eixo x (m)')
ylabel('eixo z (m)')
title('A, plano XZ (y=0)');

%Gráfico XY de A em z = 0
figure(2);
[X,Y] = meshgrid(x,y);
quiver(X, Y, squeeze(Ax(:,:,zzero))', squeeze(Ay(:,:,zzero))');
xlabel('eixo x (m)')
ylabel('eixo y (m)')
title('A, plano XY (z=0)');

%% Descobrindo valor de B

xn = (x(2:end-1));
yn = (y(2:end-1));
zn = (z(2:end-1));

Bx = zeros (length(xn),length(yn),length(zn)); %inicializa a densidade de fluxo magnético B, componente x
By = zeros (length(xn),length(yn),length(zn)); %componente y
Bz = zeros (length(xn),length(yn),length(zn)); %componente z

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
%% Grafico B

xzero2=ceil(length(xn)/2);
yzero2=ceil(length(yn)/2);
zzero2=ceil(length(zn)/2);

%Gráfico XZ de B em y = 0
figure (3)
[X,Z] = meshgrid(xn,zn);
quiver(X,Z, squeeze(Bx(:,yzero2,:))', squeeze(Bz(:,yzero2,:))')
xlabel('eixo x (m)')
ylabel('eixo z (m)')
title('B, plano XZ (y=0)');

%% Resultado no ponto

%A no ponto (0; 0; 7,3)
disp(By(xzero2,yzero2,zzero2*1.5)/u0);
