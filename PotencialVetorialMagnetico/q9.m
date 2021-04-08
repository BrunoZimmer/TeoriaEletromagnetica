% Encontre o potencial vetorial magnético de um segmento finito de fio reto 
% pelo qual passa a corrente I. Coloque o fio no eixo z, de z1 a z2. A partir 
% do potencial vetorial magnético calcule a densidade de fluxo magnético. Afim 
% de avaliar a sua reposta  calcule a magnitude do potencial vetorial magnético 
% num ponto P( x= 9,3 m, y =0, z =0), considerando I = 6,1 A, z1 = 3,2 m e z2 = 7,3 m.
clc
clear all;
close all;

%% Variáveis iniciais

u0 = 4*pi*10^-7; % permeablidade do espaço livre (em H/m)
z1 = 3.2;
z2 = 7.3;
Px = 9.3;
L = z2-z1; %comprimento do filamento (em m)
I = 6.1; %corrente no filamento

%% Espaços de calculos

%limites das coordendas x, y, z onde será calculado A
xmax = 2*Px;
ymax = L;
zmax = L;

%passos
dx = Px/10;
dy = L/20;
dz = L/20;

%posições onde será calculado A
x = -xmax:dx:xmax;
y = -ymax:dy:ymax;
z = -2*zmax:dz:2*zmax;

%índice do valor médio dos vetores x, y e z,
xzero=ceil(length(x)/2);
yzero=ceil(length(y)/2);
zzero=ceil(length(z)/2);

%vou dividir o filamento em segmentos de tamanho dzl
dzl=L/20; %passo

%posições centrais dos segmentos do filmaento
zl = z1+dzl/2:dzl:z2-dzl/2;

%% Descobrindo valor de A

%zerando o potencial vetorial magnético
Ax = zeros(length(x),length(y),length(z));
Ay = zeros(length(x),length(y),length(z));
Az = zeros(length(x),length(y),length(z));

for i = 1:length(x) % varre a coordenada x onde será calculado A
    disp(i)
    for j = 1:length(y)% varre a coordenada y onde será calculado A
        for k = 1 : length(z)% varre a coordenada z onde será calculado A
            for kl=1:length(zl) % varre a coordenada zl onde está a carga
                
                r = [x(i),y(j),z(k)]; % vetor posição onde calculamos A
                rl= [0,0,zl(kl)]; % vetor posição de um segmento dzl
                
                I_dl = I*dzl*[0,0,1]; %produto Idl, % dl = dz az
                dA = I_dl/sqrt((r-rl)*(r-rl)'); % contribuição para o potencial A devido a a corrente I no segmento dzl
                
                Ax(i,j,k)=Ax(i,j,k)+dA(1); %soma em Ax a componente x da contribuição em dA
                Ay(i,j,k)=Ay(i,j,k)+dA(2);%soma em Ay a componente y da contribuição em dA
                Az(i,j,k)=Az(i,j,k)+dA(3);%soma em Az a componente z da contribuição em dA
                
            end
        end
    end
end
Ax=(u0/(4*pi))*Ax;
Ay=(u0/(4*pi))*Ay;
Az=(u0/(4*pi))*Az;

%% Graficos

%Gráfico de A no plano y = 0
figure(1);
[X,Z] = meshgrid(x,z);
quiver(X, Z, squeeze(Ax(:,yzero,:))', squeeze(Az(:,yzero,:))');
xlabel('eixo x (m)')
ylabel('eixo z (m)')
title('A, plano XZ (y=0)');

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

xzero2=ceil(length(xn)/2);
yzero2=ceil(length(yn)/2);
zzero2=ceil(length(zn)/2);

%% Grafico

%B no plano z =0
figure
[X,Y] = meshgrid(xn,yn);
quiver(X,Y,squeeze(Bx(:,:,zzero2))' , squeeze(By(:,:,zzero2))') 
xlabel('eixo x (m)')
ylabel('eixo y (m)')

%% Resultado no ponto
%A no ponto (3,4; 0; 0)

disp(Az(int64(1.5*xzero),yzero,zzero));
%RESP =  5,51441470e-7 Wb/m.