% Determine, numericamente, a indutância mútua entre um longo fio reto e 
% uma espira em formato de triânguilo equilátero de labo b (em metros).
% Um dos vértices do triângulo é o ponto mais próximo do fio, e está a uma 
% distância d (em metros), conforme a figura. Para avaliar sua resposta,
% considere b=1 cm e d=17 cm.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%definir o cabo
%definir o triangulo

% Descobrir a densidade de fluxo B no fio por Biot-Savart
    % integrando em z em todo o espaço

% Descobrir o fluxo magnético a partir do campo B do fio no triangulo
    % integrando em x e z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

%% Variaveis iniciais

u0 = 4*pi*10^(-7); %permeabilidade magnética do espaço livre (em H/m)

b = 0.27; %lado triângulo (em m)
d = 0.21; %distância (em m)

I = 1; %corrente arbitrária (em A)

D = 2*d; %limite para as coordenadas x e y
passo = d/10;

%% Definindo espaço
%variação das coordenadas onde será calculado B
dx = passo;
dy = passo;
dz = passo;

x = -D:dx:D;
y = -D:dy:D;
z = -d/2:dz:d/2;

%variação das coordenadas do fio
xl1 = 0;
yl1 = 0;
zl1 = -2*D:dz:2*D;

%variação das coordenadas triângulo
h = b*sqrt(3)/2; %altura do triângulo
xl2 = d:dx:d+h;
yl2 = 0;
zl2 = -b/2:dz:b/2;

%índices para x = y = z = 0
xmedio = ceil(length(x)/2);
ymedio = ceil(length(y)/2);
zmedio = ceil(length(z)/2);

%% Calcular a densidade de fluxo magnético B no fio

%inicializa a densidade de fluxo magnético B
B(:,:,:,:) = zeros(3,length(x),length(y),length(z));

for i = 1:length(x) %varre as coordenadas onde B será calculado
    disp(i)
    for j = 1:length(y)
        for k = 1:length(z)
            
            for n = 1:length(zl1) %varre a coordenada do fio

                r = [x(i),y(j),z(k)];  %vetor posição do ponto de campo (onde calculamos B)
                rl1 = [xl1,yl1,zl1(n)]; %vetor posição do fio

                if ((r-rl1)*(r-rl1)' > dz/2) %evitando a divisão por zero, no caso em que r = rl

                    dL = [0,0,dz];
                    dL_x_rrl = cross(dL, r-rl1);
                    dB_ijk = I*(dL_x_rrl') / (sqrt((r-rl1)*(r-rl1)')^3); %essa é a contribuição dB devida a I

                    B(1,i,j,k) = B(1,i,j,k) + dB_ijk(1);
                    B(2,i,j,k) = B(2,i,j,k) + dB_ijk(2);
                    B(3,i,j,k) = B(3,i,j,k) + dB_ijk(3);
                end
            end
        end
    end
end
B = B*u0/(4*pi);

%% calcular o fluxo magnético do fio no triangulo

psi = 0; 
for i = 1:length(x) %varre as coordenadas onde Psi será calculado
    disp(i)
    for k = 1:length(z)

        if( (z(k) <= 1/sqrt(3)*(x(i)-d)) && (z(k) >= -1/sqrt(3)*(x(i)-d)) && (x(i) >= d) && (x(i)<= d+h) )
            
            dS = [0,dx*dz,0];
            psi = psi + dot(B(:,i,ymedio,k),dS);
            
        end
    end
end

L = psi/I;
disp(L);

