% Uma l�mina de corrente K flui na regi�o  -a < y < a no plano z =0. Calcule 
% numericamente o potencial vetorial magn�tico A e a partir desse a densidade 
% de fluxo magn�tico B, em qualquer posi��o do espa�o. Fa��o uma representa��o 
% gr�fica de A e de B no plano xz
% Afim de avaliar sua resposta calcule a componente y da intensidade de campo 
% magn�tico na posi��o (x = 0, y = 0, z = 7,3 m), considere K = 7,5 ax A/m e 
% a = 3 m.

clc
clear all
close all

%% Vari�veis iniciais

u0 = 4*pi*10^-7; % permeablidade do espa�o livre (em H/m)
a = 3; % limite da l�mina
k1 = 7.5; %na dire��o ax
zmax = 7.3;

%% Espa�os de calculos

D=2*a; % limite para as coordenadas x, y e z
passo=zmax/20; % passo

dx=passo;
dy=passo;
dz=passo;

x= -D:dx:D; %variacao da coordenada x 
y= -D:dy:D; %variacao da coordenada y
z= -2*zmax:dz:2*zmax; %variacao da coordenada z

xl= -2*D:dx:2*D; %variacao da coordenada x onde est� a carga
yl= -a:dy:a; %variacao da coordenada y onde est� a carga

%% Descobrindo valor de A

%zerando o potencial vetorial magn�tico
Ax = zeros(length(x),length(y),length(z));
Ay = zeros(length(x),length(y),length(z));
Az = zeros(length(x),length(y),length(z));

for i = 1:length(x)% varre a coordenada x onde H ser� calculado
    disp(i)
    for j = 1: length(y)  % varre a coordenada y onde H ser� calculado
        for k = 1:length(z) % varre a coordenada z onde H ser� calculado
            
            for m = 1:length(xl)  % varre a coordenada xl, onde existe a corrente
                for n = 1:length(yl)
                    
                    %vetor posi��o do ponto no espa�o
                    r = [x(i),y(j),z(k)];  
                    %vetor posi��o do ponto na carga
                    rl = [xl(m), yl(n), 0];
                    
                    K = k1;%ax
                    dS = dx*dy;
                    K_dS = K*dS*[1,0,0];
                    
                    % pequena contribuicao para campo A de K em cada posi��o
                    dA = K_dS/sqrt( (r-rl)*(r-rl)');
                    
                    Ax(i,j,k)=Ax(i,j,k)+dA(1); %soma em Ax a componente x da contribui��o em dA
                    Ay(i,j,k)=Ay(i,j,k)+dA(2);%soma em Ay a componente y da contribui��o em dA
                    Az(i,j,k)=Az(i,j,k)+dA(3);%soma em Az a componente z da contribui��o em dA
                
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

%Gr�fico XZ de A em y = 0
figure(1);
[X,Z] = meshgrid(x,z);
quiver(X, Z, squeeze(Ax(:,yzero,:))', squeeze(Az(:,yzero,:))');
xlabel('eixo x (m)')
ylabel('eixo z (m)')
title('A, plano XZ (y=0)');

%Gr�fico XY de A em z = 0
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

Bx = zeros (length(xn),length(yn),length(zn)); %inicializa a densidade de fluxo magn�tico B, componente x
By = zeros (length(xn),length(yn),length(zn)); %componente y
Bz = zeros (length(xn),length(yn),length(zn)); %componente z

%B a partir do rotacional de A 
%rot A = (dAz/dy-dAy/dz) ax + (dAx/dz-dAz/dx) ay + (dAy/dx-dAx/dy) az

for i = 2:length(x)-1% varre a coordenada x onde B ser� calculado
    disp(i)
    for j = 2 :length(y)-1 % varre a coordenada y onde B ser� calculado
        for k = 2: length(z)-1  % varre a coordenada z onde E ser� calculado
            
            dAz_dy= (Az(i,j+1,k)-Az(i,j-1,k))/2/dy; %derivada de Az em rela��o a y
            dAy_dz = (Ay(i,j,k+1)-Ay(i,j,k-1))/2/dz; %derivada de Ay em rela��o a z
            
            dAx_dz=(Ax(i,j,k+1)-Ax(i,j,k-1))/2/dz; %derivada de Ax em rela��o a z
            dAz_dx=(Az(i+1,j,k)-Az(i-1,j,k))/2/dx;%derivada de Az em rela��o a x
            
            dAy_dx=(Ay(i+1,j,k)-Ay(i-1,j,k))/2/dx; %derivada de Ay em rela��o a x
            dAx_dy=(Ax(i,j+1,k)-Ax(i,j-1,k))/2/dy;%derivada de Ax em rela��o a y
            
            Bx(i-1,j-1,k-1) = dAz_dy-dAy_dz; %a componente x de B � dAz/dy-dAy/dz
            By(i-1,j-1,k-1) = dAx_dz-dAz_dx;%a componente y de B � dAx/dz-dAz/dx
            Bz(i-1,j-1,k-1) = dAy_dx-dAx_dy;%a componente z de B � dAy/dx-dAx/dy
        end
    end
end
%% Grafico B

xzero2=ceil(length(xn)/2);
yzero2=ceil(length(yn)/2);
zzero2=ceil(length(zn)/2);

%Gr�fico XZ de B em y = 0
figure (3)
[X,Z] = meshgrid(xn,zn);
quiver(X,Z, squeeze(Bx(:,yzero2,:))', squeeze(Bz(:,yzero2,:))')
xlabel('eixo x (m)')
ylabel('eixo z (m)')
title('B, plano XZ (y=0)');

%% Resultado no ponto

%A no ponto (0; 0; 7,3)
disp(By(xzero2,yzero2,zzero2*1.5)/u0);
