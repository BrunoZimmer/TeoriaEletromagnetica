% Calcule, numericamente, a auto-indut�ncia de uma espira de raio a.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%definir uma espira no espaco de raio a 

%encontra o campo B dessa espira por Biot Savart integrando em phi

%encontrar o fluxo magn�tico causado pelo pr�prio campo da espira
%integrando em x e y na area da espira

%definir a auto indutancia por L = psi/I

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

%% Vari�veis iniciais
u0=4*pi*10^(-7); % Permeabilidade magnética do espaço livre (em H/m)

a = 1; % raio do anel

I=1; % corrente que flui no anel

%% Espa�o

D=2*a; % limite para as coordenadas x, y e z
passo=a/10; % passo

dx=passo;
dy=passo;
dz=passo;

x= -D:dx:D; %variação da coordenada x onde será calculado B
y= -D:dy:D; %variação da coordenada y onde será calculado B
z= -D:dz:D; %variação da coordenada z onde será calculado B

xmedio = ceil(length(x)/2);
ymedio = ceil(length(y)/2);
zmedio = ceil(length(z)/2);

dphi = 2*pi/14;
phil = dphi/2:dphi:2*pi-dphi/2; % notar que excluímos 2*pi, pois já temos 0.
rhol = a;

%% Calcula B

% Inicializa a densidade de fluxo magnético
B(:,:,:,:) = zeros(3,length(x),length(y),length(z));

for i = 1:length(x)% varre a coordenada x onde B será calculado
    disp(i)
    for j = 1: length(y)  % varre a coordenada y onde B será calculado
        for k = 1:length(z) % varre a coordenada z onde B será calculado
            
            for m = 1:length(phil)  % varre a coordenada phil, onde existe a corrente
                
                r = [x(i),y(j),z(k)];  %vetor posição do ponto de campo (onde calculamos B)
                
                [rl(1), rl(2), rl(3)] = pol2cart(phil(m), a, 0); % convertemos rl para coord retangulares
                
                [dL(1), dL(2), dL(3)] = pol2cart(pi/2+phil(m), a*dphi, 0);
                
                % Agora calculamos o produto vetorial dL x (r-rl)
                
                dl_x_rrl = cross(dL, r-rl);
                              
                % Utilizamos a condição abaixo para evitar a divisão por zero, no caso em que r e rl são iguais.
                if ((r-rl)*(r-rl)' > (dz/2)^2)
                    dB_ijk = dl_x_rrl'/(sqrt((r-rl)*(r-rl)')^3); % essa é a contribuição dB devida a esse dL
                    
                    B(1,i,j,k) = B(1,i,j,k) + dB_ijk(1);
                    B(2,i,j,k) = B(2,i,j,k) + dB_ijk(2);
                    B(3,i,j,k) = B(3,i,j,k) + dB_ijk(3);
                end
            end
        end
    end
end
B = (B*u0*I)/(4*pi);

%% Plots B

figure(1);
[X,Z] = meshgrid(x,z);
quiver(X, Z, squeeze(B(1,:,ymedio,:))', squeeze(B(3,:,ymedio,:))');
xlabel('eixo x (m)')
ylabel('eixo z (m)')
title('B, plano XZ (y=0)');

%% Calcular fluxo magn�tico

psi = 0;
for i = 1:length(x) %varre as coordenadas onde Psi ser� calculado
    disp(i)
    for j = 1:length(y)
        
        dS = [0,0,dx*dy];
        
        if((sqrt(x(i)^2+y(j)^2) <= a) ) %se est� na espira 1
            
            psi = psi + dot(B(:,i,j,zmedio),dS);
        end
    end
end

L = psi/I;
disp(L)