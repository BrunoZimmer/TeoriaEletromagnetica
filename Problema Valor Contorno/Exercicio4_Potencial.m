% Uma esfera de raio a contém uma densidade volumétrica de carga uniforme ?0.
% Encontre, numericamente, a energia total armazenada de duas formas (usando 
% a densidade de carga e o potencial e usando o campo elétrico devido a 
% distribuição de carga). Avalie sua resposta considerando  ?0 = 10 ?C/m³ e 
% a = 8,5 m.

%CALCULAR POTENCIAL 
clc
clear all
close all

%% Variaveis inicias
e0=8.854*10^-12; %permissividade do espaço livre (em F/m)
raio= 8.2; %raio

pv = 2.7*10^-6; %densidade volumetrica de carga (em C/m^3)

passe=raio/6;
lim=2*raio;

%% Limites de campos e cargas

x= -lim:passe:lim; %variação da coordenada x onde será calculado V
y= -lim:passe:lim; %variação da coordenada y onde será calculado V
z= -lim:passe:lim; %variação da coordenada z onde será calculado V

meio = int64(length(x)/2);%origem

dsph = pi/18;  %passe "geral"
dphi = dsph;   %passe  phi
dr= raio/8;    %passe  raio
dtheta = dsph; %passe  theta

rrl= dr/2:dr:raio-dr/2; % variação da coordenada r carga
thetal= dtheta/2:dtheta:pi-dtheta/2; % variação da coordenada theta carga
phil= dphi/2:dphi:2*pi-dphi/2; % variação da coordenada phi carga

%% Potencial a partir de pv no espaço
V(:,:,:) = zeros (length(x), length(y), length(z));

for i = 1:length(x)% varre a coordenada x onde v será calculado
    disp(i);
    for j = 1: length(y)  % varre a coordenada z onde V será calculado
        for k = 1:length(z) % varre a coordenada z onde V será calculado
            
            for m = 1:length(rrl)  % varre a coordenada rl da carga
                for n=1:length(phil) % varre a coordenada  phil da carga
                    for o=1:length(thetal) %varre a coordenada  thetal da carga
                        r = [x(i),y(j),z(k)]; %vetor posição apontando para onde estamos calculando V
                        
                        %[x,y,z] = sph2cart(azimuth,elevation,r)
                        [rl(1),rl(2),rl(3)]= sph2cart( thetal(o),phil(n),rrl(m));% vetor posição apontando para um elemento de volume dV.
                        
%                         dV = rrl(m)^2*sin(thetal(o))*dphi*pi/180*dr*dtheta*pi/180;
                        dV = rrl(m)^2*sin(thetal(o))*dphi*dr*dtheta; 
                        %dV = rrl(m)^2*sin(thetal(o))*dll*dll*dll; 
                        
                        if (sqrt((r-rl)*(r-rl)')>0.01)
                            V(i,j,k) = V(i,j,k)  + (1/(4*pi*e0))*((pv*dV)/sqrt((r-rl)*(r-rl)'))';
                        end
                        
                    end
                end
            end
        end
    end
end
%% Plot gráficos
%Grafico em Y x Z
[Y,Z] = meshgrid(y,z);
figure
[C,h] =  contour(Y,Z,squeeze(V(meio,:,:)),20); 
set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
xlabel('eixo y (m)')
ylabel('eixo z (m)')

%Grafico em X x Z
[X,Z] = meshgrid(x,z);
figure
contour(X,Z,squeeze(V(:,meio,:)),20); 
set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
xlabel('eixo z (m)')
ylabel('eixo x (m)')

%% Energia total intv(pvV)
We=0; %energia total armazenada

for i = 1:length(x)% varre a coordenada x onde v será calculado
    disp(i);
    for j = 1: length(y)  % varre a coordenada z onde V será calculado
        for k = 1:length(z) % varre a coordenada z onde V será calculado
            
            if(sqrt(x(i)^2+y(j)^2+z(k)^2)<raio) %limite da esfera 
                We = We + (0.5)*pv *V(i,j,k)*(passe^3);%formula de energia total
            end
        end
    end
end

disp(We); % resposta 419,83 kJ  - 8,5m e 10uC/m³