%enunciado

%CALCULAR POTENCIAL 
clc
clear all
close all

%% Variaveis inicias
e0 = 8.854*10^-12; %permissividade do espaço livre (em F/m)
e=1 * e0;

a = 1; %raio a que não tem valor então eu atribui 1m
L = 1;
P0 = 1;
psb = P0; % Porque o psb é igual a P na direção normal e não temos esse valor

passe = a/10;
lim = 3;

%% Limites de campos e cargas

x= -lim:passe:lim; %variação da coordenada x onde será calculado V
y= -lim:passe:lim; %variação da coordenada y onde será calculado V
z= -lim:passe:lim; %variação da coordenada z onde será calculado V

meio = int64(length(x)/2);%origem

dtheta = pi/18; %passe  theta
dp = passe; %passe  rho
%z vai ser constante porque só vou calcular as duas cascas


thetal= dtheta/2:dtheta:2*pi-dtheta/2; % variação da coordenada theta carga
pl = dp/2:dp:a;

%% Potencial a partir de pv no espaço
V(:,:,:) = zeros (length(x), length(y), length(z));

for i = 1:length(x)% varre a coordenada x onde v será calculado
    disp(i);
    for j = 1: length(y)  % varre a coordenada z onde V será calculado
        for k = 1:length(z) % varre a coordenada z onde V será calculado
            
            for n=1:length(pl) % varre a coordenada  phil da carga
                for o=1:length(thetal) %varre a coordenada  thetal da carga
                    %fazer 2 vezes a integral das tampas
                    r = [x(i),y(j),z(k)]; %vetor posição apontando para onde estamos calculando V

                    [rl(1),rl(2),rl(3)]= pol2cart( thetal(o),pl(n),L);% vetor posição apontando para um elemento de volume dV.
                    [rll(1),rll(2),rll(3)]= pol2cart( thetal(o),pl(n),0);% vetor posição apontando para um elemento de volume dV.
                    
                    dV = pl(n)*dtheta*dp; %dV em cilindricas mas com r fixo

                    if (sqrt((r-rl)*(r-rl)') > a*0.01)
                        V(i,j,k) = V(i,j,k)  + (1/(4*pi*e))*((psb*dV)/sqrt((r-rl)*(r-rl)'))';
                    end
                    
                    if (sqrt((r-rl)*(r-rl)') > a*0.01)
                        V(i,j,k) = V(i,j,k)  + (1/(4*pi*e))*((psb*dV)/sqrt((r-rll)*(r-rll)'))';
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

D = e0*E+P0;
%% Plot gráficos

% %Surface em X x Z de E
[X,Z] = meshgrid(xn,zn);
figure(1)
surf(X,Z,squeeze(E(3,:,meio,:)));
colormap parula

[X,Z] = meshgrid(xn,zn);
figure(2)
[C,h] =  contour(X,Z,squeeze(E(3,:,meio,:)),10); 
set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
xlabel('eixo z (m)')
ylabel('eixo x (m)')

hold on

quiver(X,Z,squeeze(E(1,:,meio,:))',squeeze(E(3,:,meio,:))'); 
xlabel('eixo z (m)')
ylabel('eixo x (m)')



% %Surface em X x Z de D
[X,Z] = meshgrid(xn,zn);
figure(3)
surf(X,Z,squeeze(D(3,:,meio,:)));
colormap parula

[X,Z] = meshgrid(xn,zn);
figure(4)
[C,h] =  contour(X,Z,squeeze(D(3,:,meio,:)),10); 
set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
xlabel('eixo z (m)')
ylabel('eixo x (m)')

hold on

quiver(X,Z,squeeze(D(1,:,meio,:))',squeeze(D(3,:,meio,:))'); 
xlabel('eixo z (m)')
ylabel('eixo x (m)')
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