% Uma longa fita condutora met�lica estende-se de -d/2 <y < d/2 ao longo do
% eixo z. Nessa fita, flui uma corrente I, no sentido positivo de z,
% distribu�da uniformemente. Um fio condutor se estende paralelo ao eixo z,
% passando pelo ponto y = 0, x = b. Uma corrente I passa por esse condutor 
% no sentido negativo de z. Determine numericamente a densidade de fluxo 
% magn�tico resultante em todo o espa�o. Trace o gr�fico de um plano 
% perpendicular ao comprimento do fio e da fita.
% 
% Em seguida, determine numericamente as for�as diferenciais por unidade de 
% comprimento existentes na regi�o em que existe a fita e o fio. Integre 
% numericamente essas for�as para determinar as for�as por unidade de 
% comprimento sofridas tanto pela fita quanto pelo fio.

clc
clear all
close all

%% Vari�veis iniciais

u0 = 4*pi*10^-7; % permeablidade do espa�o livre (em H/m)
d = 0.6;
b = 0.8;
I = 0.4;

%% Espa�os de calculos

passo=d/5; % passo

dx=passo;
dy=passo;
dz=passo;

x= -2*b:dx:2*b; %variacao da coordenada x
y= -d:dy:d; %variacao da coordenada y
z= -d/2:dz:d/2; %variacao da coordenada z

xmedio = ceil(length(x)/2);
ymedio = ceil(length(y)/2);
zmedio = ceil(length(z)/2);

yl1= -d/2:dy:d/2; %variacao da coordenada y onde est� a fita
zl1= -2*d:dz:2*d; %variacao da coordenada z onde est� a fita

xl2 = b;
zl2= -2*d:dy:2*d; %variacao da coordenada x onde est� a fio

%% Descobrindo valor de B

% Determine a densidade de fluxo magn�tico

%zerando a densidade de fluxo magn�tico
B(:,:,:,:) = zeros(3,length(x),length(y),length(z));
B1(:,:,:,:) = zeros(3,length(x),length(y),length(z));
B2(:,:,:,:) = zeros(3,length(x),length(y),length(z));

K1e(:,:,:,:) = zeros(3,length(x),length(y),length(z));
K2e(:,:,:,:) = zeros(3,length(x),length(y),length(z));

for i = 1:length(x)% varre a coordenada x onde H ser� calculado
    disp(i)
    for j = 1: length(y)  % varre a coordenada y onde H ser� calculado
        for k = 1:length(z) % varre a coordenada z onde H ser� calculado
            
            if (y(j) >= (-d/2) && y(j) <= (d/2) && x(i) > -dy/2 && x(i) < dy/2)
                K1e(3,i,j,k) =  I/d;
            end
            for m = 1:length(yl1)  % varre a coordenada yl, onde existe a corrente
                for n = 1:length(zl1)
                    
                    %vetor posi��o do ponto no espa�o
                    r1 = [x(i),y(j),z(k)];  
                    %vetor posi��o do ponto na carga
                    rl1 = [0, yl1(m), zl1(n)];
                    
                    K1 = [0,0,I/d]; 
                    dS1 = dy*dz; %componente de espa�o bidimensional cart.
                    
                    %produto evetorial d K x ar-rl
                    K1_x_rrl = cross(K1, r1-rl1);
                    
                    if  ((r1-rl1)*(r1-rl1)' > dz/2)     
                      % pequena contribuicao para campo A de K em cada posi��o
                      dB1_ijk = K1_x_rrl'/(sqrt((r1-rl1)*(r1-rl1)')^3); 
                
                      B1(1,i,j,k) = B1(1,i,j,k) + dB1_ijk(1)*dS1; 
                      B1(2,i,j,k) = B1(2,i,j,k) + dB1_ijk(2)*dS1;
                      B1(3,i,j,k) = B1(3,i,j,k) + dB1_ijk(3)*dS1;
                      
                    end
                end
            end
            
            if (y(j) >= (-dy/2) && y(j) <= (dy/2) && x(i) > -dx/2 && x(i) < dx/2)
                K1e(3,i,j,k) =  -I;
            end
            for n = 1:length(zl2)

                %vetor posi��o do ponto no espa�o
                r2 = [x(i),y(j),z(k)];  
                %vetor posi��o do ponto na carga
                rl2 = [b, 0, zl2(n)];

                IdL2 = [0,0,-I]; %densidade de corrente superficial em az
                dS2 = dz; %componente de espa�o bidimensional cart.

                %produto evetorial d K x ar-rl
                IdL2_x_rrl = cross(IdL2, r2-rl2);

                if  ((r2-rl2)*(r2-rl2)' > dz/2)     
                      % pequena contribuicao para campo A de K em cada posi��o
                      dB2_ijk = IdL2_x_rrl'/(sqrt((r2-rl2)*(r2-rl2)')^3); 

                      B2(1,i,j,k) = B2(1,i,j,k) + dB2_ijk(1)*dS2; 
                      B2(2,i,j,k) = B2(2,i,j,k) + dB2_ijk(2)*dS2;
                      B2(3,i,j,k) = B2(3,i,j,k) + dB2_ijk(3)*dS2;

                end
            end
        end
    end
end

B1=(u0/(4*pi))*B1;
B2=(u0/(4*pi))*B2;
B = B1+B2;


%% Graficos A
%Gr�fico XZ de A em y = 0 // Densidade de Fluxo
figure(1);
[X,Y] = meshgrid(x,y);
quiver(X, Y, squeeze(B(1,:,:,zmedio))', squeeze(B(2,:,:,zmedio))');
hold on
contour(X, Y, squeeze(B(3,:,:,zmedio))', 20);
axis([-1 1 -1 1]);
xlabel('eixo x (m)')
ylabel('eixo y (m)')
title('B, plano XY (z=0)');

figure(2);
[X,Y] = meshgrid(x,y);
quiver(X, Y, squeeze(B1(1,:,:,zmedio))', squeeze(B1(2,:,:,zmedio))');
hold on
contour(X, Y, squeeze(B1(3,:,:,zmedio))', 20);
axis([-1 1 -1 1]);
xlabel('eixo x (m)')
ylabel('eixo y (m)')
title('B1, plano XY (z=0)');

%% Descobrindo valor de For�a diferencial por unidade de comprimento

dF(:,:,:,:) = zeros(3,length(x),length(y),length(z));
dF1(:,:,:,:) = zeros(3,length(x),length(y),length(z));
dF2(:,:,:,:) = zeros(3,length(x),length(y),length(z));

dF1 = cross(K1e, B)*dS1;
dF2 = cross(K2e, B)*dS2;

dF = dF1 + dF2;

%% gr�fico da For�a Diferencial
% de uma placa
figure (3)
[X,Y] = meshgrid(x,y);
quiver(X, Y, squeeze(dF1(1,:,:,zmedio))', squeeze(dF1(2,:,:,zmedio))');
hold on
contour(X, Y, squeeze(dF1(3,:,:,zmedio))', 20);
axis([-1 1 -1 1]);
xlabel('eixo x (m)');
ylabel('eixo y (m)');
title('dF1, plano XY (z=0)');

figure (4)
[X,Y] = meshgrid(x,y);
quiver(X, Y, squeeze(dF(1,:,:,zmedio))', squeeze(dF(2,:,:,zmedio))');
hold on
contour(X, Y, squeeze(dF(3,:,:,zmedio))', 20);
axis([-1 1 -1 1]);
xlabel('eixo x (m)')
ylabel('eixo y (m)')
title('dF, plano XY (z=0)');

figure (5)
[X,Y] = meshgrid(x,y);
quiver(X, Y, squeeze(dF2(1,:,:,zmedio))', squeeze(dF2(2,:,:,zmedio))');
hold on
contour(X, Y, squeeze(dF2(3,:,:,zmedio))', 20);
axis([-1 1 -1 1]);
xlabel('eixo x (m)');
ylabel('eixo y (m)');
title('dF2, plano XY (z=0)');
%% Resultado no ponto
% determinar a magnitude da for�a por unidade de comprimento sentida por cada uma das fitas.

F1 = zeros(3,1);
F2 = zeros(3,1);

for i = 1:length(y)% varre a coordenada x onde H ser� calculado
    disp(i)
%     for k = 1:length(z) % varre a coordenada z onde H ser� calculado

    if (y(j) >= -(d/2) && y(j) <= (d/2))
        F1(1) = F1(1) + dF1(1,xmedio,j,zmedio);
        F1(2) = F1(2) + dF1(2,xmedio,j,zmedio);
        F1(3) = F1(3) + dF1(3,xmedio,j,zmedio);
    end
 
    if ((x(i) > b-dx/2) && (x(i) < b+dx/2))
        F2(1) = F2(1) + dF2(1,i,ymedio,zmedio);
        F2(2) = F2(2) + dF2(2,i,ymedio,zmedio);
        F2(3) = F2(3) + dF2(3,i,ymedio,zmedio);
    end
 
%     end
end
        
disp(F1);  
disp(F2);      
        