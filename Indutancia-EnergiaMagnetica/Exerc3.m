% Considere um cabo coaxial composto por um condutor interno sÛlido de raio
% a separado de um condutor externo de raio b, cuja espessura È d (valores
% em metros). Calcule, numericamente, a indut‚ncia por unidade de
% comprimento dessa linha de transmiss„o,usando o mÈtodo da energia. Para
% avaliar sua resposta, considere a=5,5 cm, b = 31,0 cm e d= 3,9 cm.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%definir o cabo coaxial
% condutor interno -> 0 - a (p)
% espaÁo entre condutore a - b (p)
% condutor externo b - b+d

% d È a casca externa do cabo

% Descobrir a densidade de fluxo B por Biot-Savart

% Descobrir a energia magnetica(W) W = (1/2u)*int(B≤)dv

% Finalmente pro descobrir a indutancia L = 2W/I≤
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

%% Condicoes iniciais
u0=4*pi*10^(-7); % Permeabilidade magn√©tica do espa√ßo livre (em H/m)

a = 0.037;
b = 0.101;
d = 0.05;

I = 1; %valor suposto de corrente

D=(b+d)*1; %limite para as coordenadas x, y e z
passo=a/4; % passo

dx=passo;
dy=passo;
dz=passo;

x= -1*D:dx:1*D; %varia√ß√£o da coordenada x onde ser√° calculado B
y= -1*D:dy:1*D; %varia√ß√£o da coordenada y onde ser√° calculado B
z= -2*dz:dz:2*dz; %varia√ß√£o da coordenada z onde ser√° calculado B

xmedio=ceil(length(x)/2);
ymedio=ceil(length(y)/2);
zmedio=ceil(length(z)/2);

dphi = 2*pi/12;
dp = a/4;

phil1 = dphi/2:dphi:2*pi-dphi/2;
pl1 =dp/2:dp:a-dp/2;
z1 = -1.5*D:dz:1.5*D;

phil2 = dphi/2:dphi:2*pi-dphi/2;
pl2 = b+dp/2:dp:b+d-dp/2;
z2 = -1.5*D:dz:1.5*D;

%% Descobrir a densidade de fluxo B por Biot-Savart
% Inicializa a densidade de fluxo magn√©tico

B(:,:,:,:) = zeros(3,length(x),length(y),length(z));

for i = 1:length(x)% varre a coordenada x onde B ser√° calculado
    disp(i)
    for j = 1: length(y) % varre a coordenada y onde B sera calculado
        for k = 1:length(z) % varre a coordenada z onde B ser√° calculado
            
            r = [x(i),y(j),z(k)]; %vetor posicao do ponto de campo (onde calculamos B)
            for m = 1:length(phil1) % varre a coordenada phil, onde existe a corrente
                for n = 1:length(pl1) % varre a coordenada phil, onde existe a corrente
                    for o = 1:length(z1) % varre a coordenada phil, onde existe a corrente
                        
                        [rl1(1), rl1(2), rl1(3)] = pol2cart(phil1(m), pl1(n), z1(o)); % convertemos rl para coord retangulares
                        
                        J1 = [ 0,0,I/(pi*a^2) ]; %J = I/area da seÁ„o reta do condutor interno
                        % Agora calculamos o produto vetorial dA x (r-rl)
                        J1_x_rrl1 = cross(J1, r-rl1);
                        
                        % Utilizamos a condicao abaixo para evitar a divisao por zero, no caso em que r e rl sao iguais.
                        dV1 = pl1(n)*dphi*dp*dz;

%                             dV1 = dx*dy*dz;
                        if (((r-rl1)*(r-rl1)') > (passo/2)^2)
                            dB1_ijk = ( J1_x_rrl1' / ( sqrt((r-rl1)*(r-rl1)')^3 ) ) * dV1;
                            
                            B(1,i,j,k) = B(1,i,j,k)+ (1/(4*pi))*dB1_ijk(1);
                            B(2,i,j,k) = B(2,i,j,k)+ (1/(4*pi))*dB1_ijk(2);
                            B(3,i,j,k) = B(3,i,j,k)+ (1/(4*pi))*dB1_ijk(3);
                        end
                    end
                end
            end
            for m = 1:length(phil2) % varre a coordenada phil, onde existe a corrente
                for n = 1:length(pl2) % varre a coordenada phil, onde existe a corrente
                    for o = 1:length(z2) % varre a coordenada phil, onde existe a corrente
                        
                        [rl2(1), rl2(2), rl2(3)] = pol2cart(phil2(m), pl2(n), z2(o)); % convertemos rl para coord retangulares
                        
                        J2 = [ 0,0,-I/((pi*(b+d)^2)-(pi*b^2)) ]; %J = I/area da seÁ„o reta do condutor externo
                        % Agora calculamos o produto vetorial dA x (r-rl)
                        J2_x_rrl2 = cross(J2, r-rl2);
                        
                        % Utilizamos a condicao abaixo para evitar a divisao por zero, no caso em que r e rl sao iguais.
                        dV2 = pl2(n)*dphi*dp*dz;
%                             dV2 = dx*dy*dz;
                        if (((r-rl2)*(r-rl2)') > (passo/2)^2)
                            dB2_ijk = ( J2_x_rrl2' / ( sqrt((r-rl2)*(r-rl2)')^3 ) ) * dV2;
                            
                            B(1,i,j,k) = B(1,i,j,k)+ (1/(4*pi))*dB2_ijk(1);
                            B(2,i,j,k) = B(2,i,j,k)+ (1/(4*pi))*dB2_ijk(2);
                            B(3,i,j,k) = B(3,i,j,k)+ (1/(4*pi))*dB2_ijk(3);
                        end
                    end
                end
            end
        end
    end
end
B = B*u0;
%% Graficos B

figure(1);
[X,Y] = meshgrid(x,y);
quiver(X, Y, squeeze(B(1,:,:,zmedio))', squeeze(B(2,:,:,zmedio))');
axis([-0.18 0.18 -0.18 0.18])
xlabel('eixo x (m)')
ylabel('eixo y (m)')
title('B, plano XY (z=0)');

%gr·fico de By variando em x
figure;
plot(x, B(2,:,ymedio,zmedio)); 
xlabel('eixo x')
ylabel('By(x)')
title('By variando em x');

%% Descobrir a energia magnetica(W) W = (1/2u)*int(B≤)dv

Wm=0;
dW=zeros(length(x),length(y));
% em x e y pq queremos que seja o W por z
for i = 1:length(x)% varre a coordenada x onde B ser√° calculado
    disp(i)
    for j = 1: length(y) % varre a coordenada y onde B sera calculado
        
        r = [x(i),y(j),z(k)]; %vetor posicao do ponto de campo (onde calculamos B)
        
        %W = (1/2u)*int(B≤)dv
        
        dV = dx*dy*1;
        dW(i,j) = dot(B(:,i,j,zmedio),B(:,i,j,zmedio)) ;
        
        Wm = Wm + dW(i,j)*dV; % como B È uniforme e estamos integrandi B2 num volume posso subentender a integral em z com o 1 em dx*dy*dz e sem for para z
        %
        
        
    end
end
Wm = (1/(2*u0))*Wm;

%% Finalmente pro descobrir a indutancia L = 2W/I≤
L = ( 2*Wm ) / ( I^2 );
disp(L)