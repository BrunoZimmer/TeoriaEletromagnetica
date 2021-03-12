% Dado o campo de potencial, no espaço livre, V = 100xz/(x2 + 4) V, considerando 
% que a superfície z = 0 seja a superfície de um condutor, calcule numericamente 
% a distribuição superficial de cargas para qualquer lugar no espaço, e determine 
% a carga total na porção do condutor definida por 0 < x < 9,5 m ; -9,3 m < y < 0.

clc
clear all
close all

e0=8.854*10^-12; %permissividade do espaço livre (em F/m)
% xmax = 3.7;%-4,7741n
% ymax = 7.1;
xmax = 9.2;%-1,1107e-8
ymax = 8;
zmax = 30;

%vou dividir a placa em segmentos de dimensão dx x dy
dx = xmax/50;
dy= ymax/50;
dz = zmax/50;

x = -zmax:dx:zmax; 
y = -zmax:dy:zmax; 
z = -zmax:dz:zmax; 
%% Calcular em todo campo o potencial
V(:,:,:) = zeros (length(x), length(y), length(z));

for i = 1:length(x)
    disp(i)
    for j = 1:length(y)
        for k = 1:length(z)
            
           V(i,j,k)  = (100*x(i)*z(k))/((x(i))^2 + 4);
%             V(i,j,k) = campoV;
            
        end
    end
end

%A derivadaa tem que se interna ao vampo porque se pegar um ponto do canto
%o calculo vai ta errado
xn=(x(2:end-1));
yn=(y(2:end-1));
zn=(z(2:end-1));

%Campo eletrico todo zerado 
E(1,:,:,:) = zeros (length(xn),length(yn),length(zn)); 
E(2,:,:,:) = zeros (length(xn),length(yn),length(zn));
E(3,:,:,:) = zeros (length(xn),length(yn),length(zn));

%% Calcula E
for i = 2:length(x)-1% varre a coordenada x 
    i
    for j = 2 :length(y)-1% varre a coordenada y
        for k = 2: length(z)-1  % varre a coordenada z

            %Derivada a partir de dx/dy/dz do espaço
            E(1,i-1,j-1,k-1) = -(V(i+1,j,k)-V(i-1,j,k))/2/dx; %Calcula x
            E(2,i-1,j-1,k-1) =-(V(i,j+1,k)-V(i,j-1,k))/2/dy; % Calcula Ey
            E(3,i-1,j-1,k-1) =-(V(i,j,k+1)-V(i,j,k-1))/2/dz; %calcula Ez
        end
    end
end
   
D = E*e0;
ps = D;

%% Calculando a carga total da placa
 
carga = 0;
z0 = int64(length(zn/2));

%percorrer dentro do espaço a regiao da placa
for i = 1:length(x)  
    for j = 1:length(y)  
        %Como o nosso ps esta definido em todo espaço, faz uma condição pra
        %calcular somente no espaço da placa
        if(x(i)>0 && x(i)<xmax && y(j)>-ymax && y(j)<0)
            carga = carga + ps(3,i,j,z0)*dx*dy;
        end
    end
end

disp(carga)
