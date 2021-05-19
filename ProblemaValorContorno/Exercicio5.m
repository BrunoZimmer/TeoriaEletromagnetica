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
e0 = 8.854*10^-12; %permissividade do espaço livre (em F/m)
er = 4;
e=er * e0;

a = 0.001; %raio a
b = 0.002; %raio b

Q = 1*10^-12;

ps1 = Q/(2*pi*a*passo);% o z desse caso é só um passo(um dz) pra fins de tempo
ps2 = -Q/(2*pi*b*passo);

passe = a/15;
lim = 1.2*b;

%% Limites de campos e cargas

x= -lim:passe:lim; %variação da coordenada x onde será calculado V
y= -lim:passe:lim; %variação da coordenada y onde será calculado V
z= -lim:passe:lim; %variação da coordenada z onde será calculado V

meio = int64(length(x)/2);%origem

dsph = pi/72; %passe "geral"
% dphi = dsph; %passe  phi
dtheta = dsph; %passe  theta
dz = passe; %passe  theta

thetal= dtheta/2:dtheta:2*pi; % variação da coordenada theta carga
% phil= dphi/2:dphi:2*pi-dphi/2; % variação da coordenada phi carga
zl = 0:dz:dz;

%% Potencial a partir de pv no espaço
V(:,:,:) = zeros (length(x), length(y), length(z));

for i = 1:length(x)% varre a coordenada x onde v será calculado
    disp(i);
    for j = 1: length(y)  % varre a coordenada z onde V será calculado
        for k = 1:length(z) % varre a coordenada z onde V será calculado
            
            for n=1:length(zl) % varre a coordenada  phil da carga
                for o=1:length(thetal) %varre a coordenada  thetal da carga
                    
                    r = [x(i),y(j),z(k)]; %vetor posição apontando para onde estamos calculando V

                    [rl(1),rl(2),rl(3)]= pol2cart( thetal(o),a,zl(n));% vetor posição apontando para um elemento de volume dV.
                    [rll(1),rll(2),rll(3)]= pol2cart( thetal(o),b,zl(n));% vetor posição apontando para um elemento de volume dV.

                    dV1 = a*dtheta*dz; %dV em cilindricas mas com r fixo
                    dV2 = b*dtheta*dz; 

                    if (sqrt((r-rl)*(r-rl)') > a*0.01)
                        V(i,j,k) = V(i,j,k)  + (1/(4*pi*e))*((ps1*dV1)/sqrt((r-rl)*(r-rl)'))';
                    end

                    if (sqrt((r-rll)*(r-rll)') > b*0.01)
                        V(i,j,k) = V(i,j,k)  + (1/(4*pi*e))*((ps2*dV2)/sqrt((r-rll)*(r-rll)'))';
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

[X,Y] = meshgrid(x,y);
figure(1)
surf(X,Y,squeeze(V(:,:,meio)));
colormap parula %cor binita

%Campo eletrico todo zerado 
E(1,:,:,:) = zeros (length(xn),length(yn),length(zn)); 
E(2,:,:,:) = zeros (length(xn),length(yn),length(zn));
E(3,:,:,:) = zeros (length(xn),length(yn),length(zn));

for i = 2:length(x)-1          % varre a coordenada x 
    disp(i);
    for j = 2:length(y)-1      % varre a coordenada y
        for k = 2:length(z)-1  % varre a coordenada z

            %Derivada a partir de dx/dy/dz do espaço
            E(1,i-1,j-1,k-1) = -(V(i+1,j,k)-V(i-1,j,k))/(2*passe); %Calcula Ex
            E(2,i-1,j-1,k-1) = -(V(i,j+1,k)-V(i,j-1,k))/(2*passe); % Calcula Ey
            E(3,i-1,j-1,k-1) = -(V(i,j,k+1)-V(i,j,k-1))/(2*passe); %calcula Ez
        end
    end
end


%% Plot gráficos

% %Surface em X x Y
[X,Y] = meshgrid(xn,yn);
figure(1)
surf(X,Y,squeeze(E(3,:,:,meio)));
colormap parula

%Grafico em X x Y
[X,Y] = meshgrid(xn,yn);
figure(2)
[C,h] =  contour(X,Y,squeeze(E(3,:,:,meio)),3); 
set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
xlabel('eixo x (m)')
ylabel('eixo y (m)')

hold on

quiver(X,Y,squeeze(E(1,:,:,meio))',squeeze(E(2,:,:,meio))'); 
xlabel('eixo x (m)')
ylabel('eixo y (m)')

%% Energia total ints(psE)
We=0; %energia total armazenada

for i = 1:length(xn)% varre a coordenada x 
    disp(i);
    for j = 1:length(yn) % varre a coordenada y
        for k = 1:length(zn) % varre a coordenada z
            
            dV = passe^3;
            We = We + (0.5)*e*(sqrt((E(1,i,j,k)^2)+(E(2,i,j,k)^2)+(E(3,i,j,k)^2))^2)*dV;%formula de energia total
        end
    end
end

disp(We); 