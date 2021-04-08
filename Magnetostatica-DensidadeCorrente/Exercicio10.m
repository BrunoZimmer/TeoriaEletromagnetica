%BZUma l‚mina de corrente K flui na regi„o  -a < y < a no plano z =0. A 
%partir da densidade de corrente K, calcule numericamente a intensidade de
%campo magnÈtico gerado por essa l‚mina em qualquer posiÁ„o do espaÁo.

%Afim de avaliar sua resposta calcule a componente y da intensidade de 
%campo magnÈtico na posiÁ„o (x = 0, y = 0, z = 3,3 m), 
%considere K = 5,0 ax A/m e a = 5,8 m.BZ

clc
clear all
close all

u0=4*pi*10^(-7); % Permeabilidade magnetica do espaco livre (em H/m)

a = 5.8; 
zteste = 3.3; % tamanho em z
Kx = 5;

D=2*a; % limite para as coordenadas x, y e z
passo=a/10; % passo

dx=passo;
dy=passo;
dz=passo;

x= -2*dx:dx:2*dx; %variacao da coordenada x onde sera calculado H
y= -D:dy:D; %variacao da coordenada y onde sera° calculado H
z= -zteste:dz:zteste; %variacao da coordenada z onde sera° calculado H

% Inicializa a densidade de fluxo magnetico

Hx(:,:,:) = zeros(length(x),length(y),length(z)); 
Hy(:,:,:) = zeros(length(x),length(y),length(z));
Hz(:,:,:) = zeros(length(x),length(y),length(z));

%pontos 'zeros' do espaÁo
xmedio=ceil(length(x)/2);
ymedio=ceil(length(y)/2);
zmedio=ceil(length(z)/2);

%dimens„o da lamina de corrente
xl = -D:dx:D; %variacao da coordenada x da placa
yl = -a:dy:a; %variacao da coordenada y da placa


for i = 1:length(x)% varre a coordenada x onde B ser√° calculado
  disp(i)
  for j = 1: length(y)  % varre a coordenada y onde B ser√° calculado
    for k = 1:length(z) % varre a coordenada z onde B ser√° calculado
	
      for m = 1:length(xl)  % varre a coordenada xl, onde existe a corrente
        for n = 1:length(yl)  % varre a coordenada yl, onde existe a corrente
        
            r = [x(i),y(j),z(k)];  %vetor posicao do ponto de campo (onde calculamos B)

            rl = [xl(m), yl(n), 0];

            K = [Kx, 0, 0];

            % Agora calculamos o produto vetorial K x (r-rl)
            K_x_rrl = cross(K, r-rl);


            % Utilizamos a condicao abaixo para evitar a divisao por zero, no caso em que r e rl s√£o iguais.

            if (((r-rl)*(r-rl)') > passo/2)       
                dH_ijk = K_x_rrl'/(sqrt((r-rl)*(r-rl)')^3); % essa √© a contribui√ß√£o dB devida a esse dL

                Hx(i,j,k) = Hx(i,j,k) + dH_ijk(1)*dx*dy; 
                Hy(i,j,k) = Hy(i,j,k) + dH_ijk(2)*dx*dy;
                Hz(i,j,k) = Hz(i,j,k) + dH_ijk(3)*dx*dy;

            end
        end
      end
    end
  end
end

Hx = (Hx)/(4*pi);
Hy = (Hy)/(4*pi);
Hz = (Hz)/(4*pi);

res = Hy( xmedio, ymedio, length(z));

disp(res)%-1.67676

figure(1);

[Y,Z] = meshgrid(y,z);
quiver(Y, Z, squeeze(Hy(xmedio,:,:))', squeeze(Hz(xmedio,:,:))'); 

xlabel('eixo y (m)')
ylabel('eixo z (m)')