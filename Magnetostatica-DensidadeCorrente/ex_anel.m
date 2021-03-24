
% ENG04454 - Teoria Eletromagnética Aplicada A - 2020/2
% Prof. Giovani Bulla e Prof. Raphael Brum

% Problemas numéricos de Magnetostática
% Exemplo 1: Espira com corrente

% Considere um anel filamentar de raio a, sobre o qual flui uma corrente
% de intensidade I, no sentido anti-horário, quando visto de cima.
%. O anel está localizado no plano xy, com seu centro alinhado com
% o eixo z. Calcule, numericamente, o vetor densidade de fluxo magnético
% em qualquer posição do espaço.

clc
clear all
close all

u0=4*pi*10^(-7); % Permeabilidade magnética do espaço livre (em H/m)

a=1; % raio do anel

I=1; % corrente que flui no anel

D=2*a; % limite para as coordenadas x, y e z
passo=a/8; % passo

dx=passo;
dy=passo;
dz=passo;

x=[-D:dx:D]; %variação da coordenada x onde será calculado B
y=[-D:dy:D]; %variação da coordenada y onde será calculado B
z=[-D:dz:D]; %variação da coordenada z onde será calculado B

% Inicializa a densidade de fluxo magnético

Bx(:,:,:) = zeros(length(x),length(y),length(z)); 
By(:,:,:) = zeros(length(x),length(y),length(z));
Bz(:,:,:) = zeros(length(x),length(y),length(z));

% Índice do valor médio dos vetores de x, y e z.
% (Neste exemplo, são os índices em que x=y=z=0.)

i_xmedio=ceil(length(x)/2);
j_ymedio=ceil(length(y)/2);
k_zmedio=ceil(length(z)/2);

% A corrente existe na região de um círculo de raio a, centrado na origem,
% no plano zl=0, para qualquer ângulo phil (0 a 2*pi).

% Vamos definir os "diferenciais de corrente" I*dL em termos de phil. Como rhol
% é fixo, vamos iterar apenas em phil:

dphi=2*pi/144;
phil=[0:dphi:2*pi-dphi/2]; % notar que excluímos 2*pi, pois já temos 0.


for i = 1:length(x)% varre a coordenada x onde B será calculado
  disp(i)
  for j = 1: length(y)  % varre a coordenada y onde B será calculado
    for k = 1:length(z) % varre a coordenada z onde B será calculado
	
      for m = 1:length(phil)  % varre a coordenada phil, onde existe a corrente
        
		r = [x(i),y(j),z(k)];  %vetor posição do ponto de campo (onde calculamos B)
        
		% o vetor rl é definido em cilíndricas, com rhol=a, phil variável, z=0
		% ele aponta para onde o diferencial dL se encontra.

		[rl(1), rl(2), rl(3)] = pol2cart(phil(m), a, 0); % convertemos rl para coord retangulares
		
		% IMPORTANTE: o vetor posição do dL não deve ser confundido com si próprio. O vetor dL, diferentemente de seu
        % vetor posição, aponta tangencialmente ao círculo. No problema, ele deve apontar na direção â_phi(phil).
		% â_phi(phil) e â_rho(phil) estão defasados de 90º.  O módulo de dL é rhol*dphi = a*dphi. Assim, determinamos
		% dL e convertemos para retangulares:
		
		[dL(1), dL(2), dL(3)] = pol2cart(pi/2+phil(m), a*dphi, 0);
        
		% Agora calculamos o produto vetorial dL x (r-rl)
			
		dl_x_rrl = cross(dL, r-rl);
		      
		
		% Utilizamos a condição abaixo para evitar a divisão por zero, no caso em que r e rl são iguais.
		
		if  ((r-rl)*(r-rl)' > dz/2)       
          dB_ijk = dl_x_rrl'/(sqrt((r-rl)*(r-rl)')^3); % essa é a contribuição dB devida a esse dL
		
		  Bx(i,j,k) = Bx(i,j,k) + dB_ijk(1); 
	      By(i,j,k) = By(i,j,k) + dB_ijk(2);
  		  Bz(i,j,k) = Bz(i,j,k) + dB_ijk(3);

        end

      end
    end
  end
end

Bx = Bx*u0*I/4*pi;
By = By*u0*I/4*pi;
Bz = Bz*u0*I/4*pi;

figure(1);

[X,Z] = meshgrid(x,z);
quiver(X, Z, squeeze(Bx(:,j_ymedio,:))', squeeze(Bz(:,j_ymedio,:))'); 

xlabel('eixo x (m)')
ylabel('eixo z (m)')
title('B, plano XZ (y=0)');

figure(2);

[Y,Z] = meshgrid(y,z);
quiver(Y, Z, squeeze(By(i_xmedio,:,:))', squeeze(Bz(i_xmedio,:,:))');

xlabel('eixo y (m)')
ylabel('eixo z (m)')
title('B, plano YZ (x=0)');


