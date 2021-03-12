% Uma esfera de raio a contém uma densidade volumétrica de carga uniforme ?0.
% Encontre, numericamente, a energia total armazenada de duas formas (usando 
% a densidade de carga e o potencial e usando o campo elétrico devido a 
% distribuição de carga). Avalie sua resposta considerando  p0 = 2.7 C/m³ e 
% a = 8,2 m.

%CALCULAR POTENCIAL 
clc
clear all
close all

%% Variaveis inicias
e0=8.854*10^-12; %permissividade do espaço livre (em F/m)
raio= 8.2; %raio

pv = 2.7*10^-6; %densidade volumetrica de carga (em C/m^3)

passe=raio/5;
lim=4*raio;

%% Limites de campos e cargas

x= -lim:passe:lim; %variação da coordenada x onde será calculado V
y= -lim:passe:lim; %variação da coordenada y onde será calculado V
z= -lim:passe:lim; %variação da coordenada z onde será calculado V

meio = int64(length(x)/2);%origem

dsph = pi/18; %passe "geral"
dphi = dsph; %passe  phi
dr= raio/6; %passe  phi
dtheta = dsph; %passe  phi

rrl= dr/2:dr:raio-dr/2; % variação da coordenada r carga
thetal= dtheta/2:dtheta:pi-dtheta/2; % variação da coordenada theta carga
phil= dphi/2:dphi:2*pi-dphi/2; % variação da coordenada phi carga

%% Potencial a partir de pv no espaço
V(:,:,:) = zeros (length(x), length(y), length(z));

for i = 1:length(x)% varre a coordenada x onde v será calculado
    disp(i);
    for j = 1: length(y)  % varre a coordenada z onde V será calculado
        for k = 1:length(z) % varre a coordenada z onde V será calculado
            
            for n=1:length(phil) % varre a coordenada  phil da carga
                for o=1:length(thetal) %varre a coordenada  thetal da carga
                    for m=1:length(rrl) %varre a coordenada  thetal da carga
                        r = [x(i),y(j),z(k)]; %vetor posição apontando para onde estamos calculando V

                        [rl(1),rl(2),rl(3)]= sph2cart( thetal(o),phil(n),rrl(m));% vetor posição apontando para um elemento de volume dV.

                        dV = rrl(m)^2*sin(thetal(o))*dphi*dr*dtheta; 

                        if (sqrt((r-rl)*(r-rl)')>0.01)
                            V(i,j,k) = V(i,j,k)  + (1/(4*pi*e0))*((pv*dV)/sqrt((r-rl)*(r-rl)'))';
                        end

                    end
                end
            end
        end
    end
end

%% Campo elétrico derivando 
xn=(x(2:end-1));
yn=(y(2:end-1));
zn=(z(2:end-1));

%Campo eletrico todo zerado 
E(1,:,:,:) = zeros (length(xn),length(yn),length(zn)); 
E(2,:,:,:) = zeros (length(xn),length(yn),length(zn));
E(3,:,:,:) = zeros (length(xn),length(yn),length(zn));

for i = 2:length(x)-1% varre a coordenada x 
    disp(i);
    for j = 2 :length(y)-1% varre a coordenada y
        for k = 2: length(z)-1  % varre a coordenada z

            %Derivada a partir de dx/dy/dz do espaço
            E(1,i-1,j-1,k-1) = -(V(i+1,j,k)-V(i-1,j,k))/(2*passe); %Calcula Ex
            E(2,i-1,j-1,k-1) = -(V(i,j+1,k)-V(i,j-1,k))/(2*passe); % Calcula Ey
            E(3,i-1,j-1,k-1) = -(V(i,j,k+1)-V(i,j,k-1))/(2*passe); %calcula Ez
        end
    end
end

%% Energia total ints(psE)
We=0; %energia total armazenada

for i = 1:length(xn)% varre a coordenada x 
    disp(i);
    for j = 1:length(yn) % varre a coordenada y
        for k = 1:length(zn) % varre a coordenada z
            dV = passe^3;

            We = We + (0.5)*e0*(sqrt((E(1,i,j,k)^2)+(E(2,i,j,k)^2)+(E(3,i,j,k)^2))^2)*dV;%formula de energia total
            
        end
    end
end

disp(We); % resposta 25,319 kJ  