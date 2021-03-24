%Um disco de raio a pertence ao plano xy, com o eixo z passando pelo seu 
%centro. Uma carga superficial de densidade uniforme ps est� presente no 
%disco, que gira em volta do eixo z em uma velocidade angular de w rad/s. 
%A partir da Lei de Biot-Savart e da densidade superficial de cargas, 
%determine numericamente a intensidade de campo magn�tico em todo o espa�o.
%Para fins de confer�ncia, avalie sua resposta considerando a = 6,9 m,
%ps = 5,5 C/m2, z = 2,0 m e w = 6,4 rad/s

%B - densidade
%H - intensidade

clc
clear all
close all
%% Condicoes iniciais
u0=4*pi*10^(-7); % Permeabilidade magnética do espaço livre (em H/m)

a=6.9; % raio do disco
zteste = 2; % tamanho em z
w= 6.4; % velocidade angular em rad/s
ps = 5.5;

D=a; % limite para as coordenadas x, y e z
passo=a/15; % passo

dx=passo;
dy=passo;
dz=passo;

x= -D:dx:D; %variação da coordenada x onde será calculado B
y= -D:dy:D; %variação da coordenada y onde será calculado B
z= -2*zteste:dz:2*zteste; %variação da coordenada z onde será calculado B

% Inicializa a densidade de fluxo magnético

Hx(:,:,:) = zeros(length(x),length(y),length(z)); 
Hy(:,:,:) = zeros(length(x),length(y),length(z));
Hz(:,:,:) = zeros(length(x),length(y),length(z));

xmedio=ceil(length(x)/2);
ymedio=ceil(length(y)/2);
zmedio=ceil(length(z)/2);

dphi = 2*pi/144;
dp = a/10;

phil = 0:dphi:2*pi-dphi/2; 
pl = 0:dp:a; 


for i = 1:length(x)% varre a coordenada x onde B será calculado
  disp(i)
  for j = 1: length(y)  % varre a coordenada y onde B sera calculado
    for k = 1:length(z) % varre a coordenada z onde B será calculado
	
      for m = 1:length(phil)  % varre a coordenada phil, onde existe a corrente
        for n = 1:length(pl)  % varre a coordenada phil, onde existe a corrente
        
            r = [x(i),y(j),z(k)];  %vetor posição do ponto de campo (onde calculamos B)


            [rl(1), rl(2), rl(3)] = pol2cart(phil(m), pl(n), 0); % convertemos rl para coord retangulares

            [K(1), K(2), K(3)] = pol2cart(pi/2+phil(m),pl(n)*ps*w, 0);

            % Agora calculamos o produto vetorial dA x (r-rl)

            K_x_rrl = cross(K, r-rl);


            % Utilizamos a condição abaixo para evitar a divisão por zero, no caso em que r e rl são iguais.

            if (((r-rl)*(r-rl)') > dz/2)       
                dH_ijk = K_x_rrl'/(sqrt((r-rl)*(r-rl)')^3); 

                Hx(i,j,k) = Hx(i,j,k) + (1/(4*pi))*dH_ijk(1)*pl(n)*dphi*dp; %tirar os d
                Hy(i,j,k) = Hy(i,j,k) + (1/(4*pi))*dH_ijk(2)*pl(n)*dphi*dp;
                Hz(i,j,k) = Hz(i,j,k) + (1/(4*pi))*dH_ijk(3)*pl(n)*dphi*dp;

            end
        end
      end
    end
  end
end


res = 0;
res = res + Hx(xmedio,ymedio,ceil(zmedio*1.5));
res = res + Hy(xmedio,ymedio,ceil(zmedio*1.5));
res = res + Hz(xmedio,ymedio,ceil(zmedio*1.5));

disp(res)%65,84

% figure(1);
% 
% [X,Z] = meshgrid(x,z);
% quiver(X, Z, squeeze(Hx(:,ymedio,:))', squeeze(Hz(:,ymedio,:))');
% 
% xlabel('eixo x (m)')
% ylabel('eixo z (m)')
% title('B, plano XZ (y=0)');
% 
% figure(2);
% 
% [Y,Z] = meshgrid(y,z);
% quiver(Y, Z, squeeze(Hy(xmedio,:,:))', squeeze(Hz(xmedio,:,:))');
% 
% xlabel('eixo y (m)')
% ylabel('eixo z (m)')
% title('B, plano YZ (x=0)');
