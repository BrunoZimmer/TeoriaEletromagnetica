% Uma esfera de raio a cont�m uma densidade volum�trica de carga uniforme ?0.
% Encontre, numericamente, a energia total armazenada de duas formas (usando 
% a densidade de carga e o potencial e usando o campo el�trico devido a 
% distribui��o de carga). Avalie sua resposta considerando  ?0 = 10 ?C/m� e 
% a = 8,5 m.

%CALCULAR POTENCIAL 
clc
clear all
close all

%% Variaveis inicias
e0=8.854*10^-12; %permissividade do espa�o livre (em F/m)
raio= 6.7; %raio

ps = 8.1*10^-6; %densidade volumetrica de carga (em C/m^2)

passe=raio/10;
lim=3*raio;

%% Limites de campos e cargas

x= -lim:passe:lim; %varia��o da coordenada x onde ser� calculado V
y= -lim:passe:lim; %varia��o da coordenada y onde ser� calculado V
z= -lim:passe:lim; %varia��o da coordenada z onde ser� calculado V

meio = int64(length(x)/2);%origem

dsph = pi/18; %passe "geral"
dphi = dsph; %passe  phi
dr= raio/10; %passe  phi
dtheta = dsph; %passe  phi

rrl= raio; % varia��o da coordenada r carga
thetal= dtheta/2:dtheta:pi-dtheta/2; % varia��o da coordenada theta carga
phil= dphi/2:dphi:2*pi-dphi/2; % varia��o da coordenada phi carga

%% Potencial a partir de pv no espa�o
V(:,:,:) = zeros (length(x), length(y), length(z));

for i = 1:length(x)% varre a coordenada x onde v ser� calculado
    disp(i);
    for j = 1: length(y)  % varre a coordenada z onde V ser� calculado
        for k = 1:length(z) % varre a coordenada z onde V ser� calculado
            
%             for m = 1:length(rrl)  % varre a coordenada rl da carga
            for n=1:length(phil) % varre a coordenada  phil da carga
                for o=1:length(thetal) %varre a coordenada  thetal da carga
                    r = [x(i),y(j),z(k)]; %vetor posi��o apontando para onde estamos calculando V

                    %[x,y,z] = sph2cart(azimuth,elevation,r)
                    [rl(1),rl(2),rl(3)]= sph2cart( thetal(o),pi/2 - phil(n),rrl);% vetor posi��o apontando para um elemento de volume dV.

                    dS = rrl^2*sin(thetal(o))*dphi*dtheta; 

                    if (sqrt((r-rl)*(r-rl)')>0.01)
                        V(i,j,k) = V(i,j,k)  + (1/(4*pi*e0))*((ps*dS)/sqrt((r-rl)*(r-rl)'))';
                    end

                end
            end
%             end
        end
    end
end
%% Plot gr�ficos
Grafico em Y x Z
[Y,Z] = meshgrid(y,z);
figure
[C,h] =  contour(Y,Z,squeeze(V(meio,:,:)),20); 
set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
xlabel('eixo y (m)')
ylabel('eixo z (m)')
% 
% %Grafico em X x Z
% [X,Z] = meshgrid(x,z);
% figure
% contour(X,Z,squeeze(V(:,meio,:)),20); 
% set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
% xlabel('eixo z (m)')
% ylabel('eixo x (m)')

%% Energia total intv(pvV)
We=0; %energia total armazenada

for i = 1:length(x)% varre a coordenada x onde v ser� calculado
    disp(i);
    for j = 1: length(y)  % varre a coordenada z onde V ser� calculado
        for k = 1:length(z) % varre a coordenada z onde V ser� calculado
            dV = passe^3;
            
            if(sqrt(x(i)^2+y(j)^2+z(k)^2) < 1.04*raio) && (sqrt(x(i)^2+y(j)^2+z(k)^2) > 0.96*raio) %limite da esfera 
                We = We + (0.5)*ps *V(i,j,k)*dV;%formula de energia total
            end
        end
    end
end

disp(We); % resposta 18,413 kJ  - 8,5m e 6.5uC/m�