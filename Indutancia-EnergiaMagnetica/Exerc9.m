% Dois anÈis circulares coaxiais de raios a e b tÍm seus eixos alinhados
% com o eixo z. O anel de raio a tem centro na origem de coordenadas,
% enquanto que o anel de raio b tem seu centro em coordenadas x, y e z
% arbitr·rias.  Determine, numericamente, a indut‚ncia m˙tua entre os anÈis.
% Afim de avaliar sua resposta considere que os dois circuitos est„o no espaÁo
% livre, dispostos como mostrado na figura, com  a =1,8 mm, b = 2,1 mm e h = 26,3 mm.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%definir as espiras no espaco de raio a e de raio b

%encontra o campo B da espira de raio b por Biot Savart integrando e 
%generalizando para qualquer ponto xc,yc,zc no espaÁo

%encontrar o fluxo magnÈtico causado pelo campo da outra espira
%integrando em x e y na area da espira

%definir a auto indutancia por L = psi/I

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
%% Vari·veis fornecidas

u0=4*pi*10^(-7); % Permeabilidade magn√©tica do espa√ßo livre (em H/m)

a = 1.1*10^(-3); %raio do anel 1
b = 2.1*10^(-3); %raio do anel 2
h = 26.6*10^(-3); %dist‚ncia entre as espiras
% a = 10.0011; % raio do anel 1
% b = 20.0021; % raio do anel 2
% h = 200.0266;
% a = 0.0011; % raio do anel 1
% b = 0.0021; % raio do anel 2
% h = 0.0266;

I=1; % corrente que flui no anel

%% Vari·veis Iniciais
lim = 2*b; % limite para as coordenadas x, y e z
passo = b/10; % passo

dx = passo;
dy = passo;
dz = passo;

x = -lim:dx:lim; %varia√ß√£o da coordenada x onde ser√° calculado B
y = -lim:dy:lim; %varia√ß√£o da coordenada y onde ser√° calculado B
z = -1.5*h:dz:1.5*h; %varia√ß√£o da coordenada z onde ser√° calculado B

xmedio=ceil(length(x)/2);
ymedio=ceil(length(y)/2);
zmedio=ceil(length(z)/2);

dphi=2*pi/12;
phil=dphi/2:dphi:2*pi-dphi/2;

xc = 0;
yc = 0;
zc = h;
%% Calcula B

% Inicializa a densidade de fluxo magnetico
B(:,:,:,:) = zeros(3,length(x),length(y),length(z));

for i = 1:length(x)% varre a coordenada x onde B ser√° calculado
    disp(i)
    for j = 1: length(y)  % varre a coordenada y onde B ser√° calculado
        for k = 1:length(z) % varre a coordenada z onde B ser√° calculado
            
            for m = 1:length(phil)  % varre a coordenada phil, onde existe a corrente
                
                r = [x(i),y(j),z(k)];  %vetor posi√ß√£o do ponto de campo (onde calculamos B)
                [rl(1), rl(2), rl(3)] = pol2cart(phil(m), b, h); % convertemos rl para coord retangulares
                rl = [rl(1)+xc, rl(2)+yc, rl(3)];
                
                [dL(1), dL(2), dL(3)] = pol2cart(pi/2+phil(m), b*dphi, 0);%quero direcao phi

                dl_x_rrl = cross(dL, r-rl);
 
                if  ((r-rl)*(r-rl)' > (dz/2)^2)
                    dB = dl_x_rrl'/(sqrt((r-rl)*(r-rl)')^3); % essa √© a contribui√ß√£o dB devida a esse dL

                    B(1,i,j,k) = B(1,i,j,k) + dB(1);
                    B(2,i,j,k) = B(2,i,j,k) + dB(2);
                    B(3,i,j,k) = B(3,i,j,k) + dB(3);
                    
                end
                
            end
        end
    end
end

B = (B*u0*I)/(4*pi);

%% Calcula Fluxo

psi = 0;
for i = 1:length(x) %varre as coordenadas onde Psi ser· calculado
    disp(i)
    for j = 1:length(y)
        
        dS1 = [0,0,dx*dy];
        
        if( (sqrt(x(i)^2+y(j)^2) <= a) ) %se est· na espira 1
            %psi = integral(B*dS)           
            psi = psi + dot(B(:,i,j,zmedio),dS1);
        end
    end
end

%% Calcula L
L = psi/I;
disp(L);

%% Plot B
%gr·fico vetorial de B no plano XY
[Y,Z] = meshgrid(y,z);
figure (1);
quiver(Y, Z, squeeze(B(2,xmedio,:,:))', squeeze(B(3,xmedio,:,:))'); 
% AXIS([XMIN XMAX YMIN YMAX])
axis([-0.004 0.004 0.022 0.032])
xlabel('eixo y (m)')
ylabel('eixo z (m)')
title('B, plano YZ (x=0)');