% Um longo cilindro de raio a tem magnetização M = k?² a?, onde k é uma 
% constante. Determine numericamente o campo vetorial magnetização. Trace
% um gráfico de M em um plano paralelo às linhas de campo. A partir de M,
% determine numericamente as densidades de corrente ligadas, e trace o
% gráficos dessas densidades, escolhendo o plano de visualização mais
% conveniente. A partir dessas densidades, determine a densidade de fluxo
% magnético em todo o espaço. Trace um gráfico da densidade de fluxo
% magnético no mesmo plano escolhido para M (em um gráfico separado).

clc
clear all
close all

%% Variáveis iniciais

u0 = 4*pi*10^-7; % permeablidade do espaço livre (em H/m)
a = 0.07;
k = 6.9;
p = 0.003;

%% Espaços de calculos

passo=a/6; % passo

dx=passo;
dy=passo;
dz=passo;

x= -2*a:dx:2*a; %variacao da coordenada x
y= -2*a:dy:2*a; %variacao da coordenada y
z= -a:dz:a; %variacao da coordenada z

xmedio = ceil(length(x)/2);
ymedio = ceil(length(y)/2);
zmedio = ceil(length(z)/2);

dp = p/5;
dphi =  pi/6;

pl= dp:dp:a; %variacao da coordenada p onde está a fita
phil= 0:dphi:((2*pi)-dphi); %variacao da coordenada phi onde está a fita
zl= -3*a:dz:3*a; %variacao da coordenada z onde está a fita

%% Descobrindo valor de M

% Determine numericamente o campo vetorial magnetização - só pega e defini pra todo o espaço

%zerando a campo vetorial magnetização
M(:,:,:,:) = zeros(3,length(x),length(y),length(z));

for i = 1:length(x)% varre a coordenada x onde H será¡ calculado
    disp(i)
    for j = 1: length(y)  % varre a coordenada y onde H será calculado
        for k = 1:length(z) % varre a coordenada z onde H será calculado 
            % M = k?² a?
            
            rho = sqrt((x(i)^2)+(y(j)^2));
            
            if ( rho <= a )
                M(1,i,j,k) = M(1,i,j,k) + k*(rho^2)*(-y(j)/sqrt((x(i)^2)+(y(j)^2))); 
                M(2,i,j,k) = M(2,i,j,k) + k*(rho^2)*(x(i)/sqrt((x(i)^2)+(y(j)^2))); 
                M(3,i,j,k) = M(3,i,j,k); 
            end
            
        end
    end
end


%% Graficos M
%Gráfico XY de M em z = 0
figure(1);
[X,Y] = meshgrid(x,y);
quiver(X, Y, squeeze(M(1,:,:,zmedio))', squeeze(M(2,:,:,zmedio))');
hold on
contour(X, Y, squeeze(M(3,:,:,zmedio))', 20);
% axis([-1 1 -1 1]);
xlabel('eixo x (m)')
ylabel('eixo y (m)')
title('M, plano XY (z=0)');

%% Descobrindo valor de J

% Determine numericamente as densidades de corrente ligadas
% rotacional de M
% K = -J

%Derivada é um espaço menor tirando o "cantos" do espaço
xn = (x(2:end-1));
yn = (y(2:end-1));
zn = (z(2:end-1));

%Zerando densidades de corrente ligadas
J(:,:,:,:)  = zeros (3,length(xn),length(yn),length(zn)); 

%J a partir do rotacional de M 
%rot M = (dMz/dy-dMy/dz) ax + (dMx/dz-dMz/dx) ay + (dMy/dx-dMx/dy) az

for i = 2:length(x)-1
    disp(i)
    for j = 2 :length(y)-1 
        for k = 2: length(z)-1  
            
            if (sqrt(x(i)^2+y(j)^2) <0.96*a)
                dMz_dy= (M(3,i,j+1,k)-M(3,i,j-1,k))/2/dy; %derivada de Mz em relação a y
                dMy_dz = (M(2,i,j,k+1)-M(2,i,j,k-1))/2/dz; %derivada de My em relação a z

                dMx_dz=(M(1,i,j,k+1)-M(1,i,j,k-1))/2/dz; %derivada de Mx em relação a z
                dMz_dx=(M(3,i+1,j,k)-M(3,i-1,j,k))/2/dx;%derivada de Mz em relação a x

                dMy_dx=(M(2,i+1,j,k)-M(2,i-1,j,k))/2/dx; %derivada de My em relação a x
                dMx_dy=(M(1,i,j+1,k)-M(1,i,j-1,k))/2/dy;%derivada de Mx em relação a y

                J(1,i-1,j-1,k-1) = dMz_dy-dMy_dz; %a componente x de J é dMz/dy-dMy/dz
                J(2,i-1,j-1,k-1) = dMx_dz-dMz_dx;%a componente y de J é dMx/dz-dMz/dx
                J(3,i-1,j-1,k-1) = dMy_dx-dMx_dy;%a componente z de J é dMy/dx-dmx/dy
            end
        end
    end
end

xmedio2=ceil(length(xn)/2);
ymedio2=ceil(length(yn)/2);
zmedio2=ceil(length(zn)/2);

%% Graficos J
%Gráfico XY de M em z = 0
figure(2);
[X,Y] = meshgrid(xn,yn);
quiver(X, Y, squeeze(J(1,:,:,zmedio2))', squeeze(J(2,:,:,zmedio2))');
hold on
contour(X, Y, squeeze(J(3,:,:,zmedio2))', 5);
% axis([-1 1 -1 1]);
xlabel('eixo x (m)')
ylabel('eixo y (m)')
title('J, plano XY (z=0)');


%Gráfico XZ de M em y = 0
figure(3);
[X,Z] = meshgrid(xn,zn);
quiver(X, Z, squeeze(J(1,:,ymedio2,:))', squeeze(J(3,:,ymedio2,:))');
hold on
contour(X, Z, squeeze(J(2,:,ymedio2,:))', 20);
axis([-0.2 0.2 -0.1 0.1]);
xlabel('eixo x (m)')
ylabel('eixo z (m)')
title('J, plano XZ (y=0)');

%% Descobrindo valor de K

K = zeros (3,length(x),length(y),length(z)); %componente x

for i = 1:length(xn) % varre a coordenada x onde será calculado A
    disp(i)
    for j = 1:length(yn)% varre a coordenada y onde será calculado A
        for k = 2 : length(zn)% varre a coordenada z onde será calculado A
            
            for m=1:length(pl) % varre a coordenada pl onde está a carga
                for n=1:length(zl) % varre a coordenada zl onde está a carga
                    for o=1:length(phil) % varre a coordenada phil onde está a carga
                        
                        rhoK = sqrt(x(i)^2 + y(j)^2);
                        
                        if (rhoK >0.96*a && rhoK <1.04*a)
                            M_x_an = cross([0,k*(pl(n))^2,0],[1,0,0]);

                            K(1,i-1,j-1,k-1) = M_x_an(1); 
                            K(2,i-1,j-1,k-1) = M_x_an(2);
                            K(3,i-1,j-1,k-1) = M_x_an(3);
                        end
                    end
                end
            end
        end
    end
end

%% Graficos K
%Gráfico XY de K em z = 0
figure(4);
[X,Y] = meshgrid(x,y);
quiver(X, Y, squeeze(K(1,:,:,zmedio))', squeeze(K(2,:,:,zmedio))');
hold on
contour(X, Y, squeeze(K(3,:,:,zmedio))', 5);
% axis([-1 1 -1 1]);
xlabel('eixo x (m)')
ylabel('eixo y (m)')
title('K, plano XY (z=0)');


%Gráfico XZ de K em y = 0
figure(5);
[X,Z] = meshgrid(x,z);
quiver(X, Z, squeeze(K(1,:,ymedio,:))', squeeze(K(3,:,ymedio,:))');
hold on
contour(X, Z, squeeze(K(2,:,ymedio,:))', 20);
% axis([-1 1 -1 1]);
xlabel('eixo x (m)')
ylabel('eixo z (m)')
title('K, plano XZ (y=0)');

%% Determine o Potencial vetorial magnetico
% Determine a densidade de fluxo magnético em todo o espaço

% encontrar A a partir de J e K
%zerando o potencial vetorial magnético
A = zeros(3,length(xn),length(yn),length(zn));

for i = 1:length(xn) % varre a coordenada x onde será calculado A
    disp(i)
    for j = 1:length(yn)% varre a coordenada y onde será calculado A
        for k = 1 : length(zn)% varre a coordenada z onde será calculado A
            
            for m=1:length(pl) % varre a coordenada pl onde está a carga
                for n=1:length(zl) % varre a coordenada zl onde está a carga
                    for o=1:length(phil) % varre a coordenada phil onde está a carga

                        r = [xn(i),yn(j),zn(k)]; % vetor posição onde calculamos A

                        % rl aponta para o diferencial dL
                        [rl(1), rl(2), rl(3)] = pol2cart(phil(o), pl(m), zl(n));

                        rhoA = sqrt(x(i)^2 + y(j)^2);
                        dV = pl(m)*dp*dz*dphi;
                        
                        if(rhoA < a)
                            dA1 = (J*dV)/sqrt((r-rl)*(r-rl)'); % contribuição para o potencial A devido a a corrente I no segmento dzl

                            A(1,i,j,k)= A(1,i,j,k) + dA1(1,i,j,k); %soma em Ax a componente x da contribuição em dA
                            A(2,i,j,k)= A(2,i,j,k) + dA1(2,i,j,k);%soma em Ay a componente y da contribuição em dA
                            A(3,i,j,k)= A(3,i,j,k) + dA1(3,i,j,k);%soma em Az a componente z da contribuição em dA
                        end
                        
                        dS = pl(m)*dp*dz;
                        if (rhoK == a)%(rhoK >0.96*raio && rhoK <1.04*raio)
                            dA2 = (K*dS)/sqrt((r-rl)*(r-rl)'); % contribuição para o potencial A devido a a corrente I no segmento dzl

                            A(1,i,j,k)= A(1,i,j,k) + dA2(1,i,j,k); %soma em Ax a componente x da contribuição em dA
                            A(2,i,j,k)= A(2,i,j,k) + dA2(2,i,j,k);%soma em Ay a componente y da contribuição em dA
                            A(3,i,j,k)= A(3,i,j,k) + dA2(3,i,j,k);%soma em Az a componente z da contribuição em dA

                        end
                    end
                end
            end
        end
    end
end
A=(u0/(4*pi))*A;


%% Graficos A
%Gráfico XY de A em z = 0
figure(6);
[X,Y] = meshgrid(xn,yn);
quiver(X, Y, squeeze(A(1,:,:,zmedio))', squeeze(A(2,:,:,zmedio))');
hold on
contour(X, Y, squeeze(A(3,:,:,zmedio))', 5);
% axis([-1 1 -1 1]);
xlabel('eixo x (m)')
ylabel('eixo y (m)')
title('A, plano XY (z=0)');


%Gráfico XZ de K em y = 0
figure(7);
[X,Z] = meshgrid(x,z);
quiver(X, Z, squeeze(K(1,:,ymedio,:))', squeeze(K(3,:,ymedio,:))');
hold on
contour(X, Z, squeeze(K(2,:,ymedio,:))', 20);
% axis([-1 1 -1 1]);
xlabel('eixo x (m)')
ylabel('eixo z (m)')
title('A, plano XZ (y=0)');
%% B a partir de rotacional do A

%Derivada é um espaço menor tirando o "cantos" do espaço
xn2 = (xn(2:end-1));
yn2 = (yn(2:end-1));
zn2 = (zn(2:end-1));

%Zerando densidades de fluxo magnetico
B(:,:,:,:)  = zeros (3,length(xn2),length(yn2),length(zn2)); 

%B a partir do rotacional de A 

for i = 2:length(xn)-1
    disp(i)
    for j = 2:length(yn)-1
        disp(i)
        for k = 2: length(zn)-1  
            
            dAz_dy= (A(3,i,j+1,k)-A(3,i,j-1,k))/2/dy; %derivada de Az em relação a y
            dAy_dz = (A(2,i,j,k+1)-A(2,i,j,k-1))/2/dz; %derivada de Ay em relação a z
            
            dAx_dz=(A(1,i,j,k+1)-A(1,i,j,k-1))/2/dz; %derivada de Ax em relação a z
            dAz_dx=(A(3,i+1,j,k)-A(3,i-1,j,k))/2/dx;%derivada de Az em relação a x
            
            dAy_dx=(A(2,i+1,j,k)-A(2,i-1,j,k))/2/dx; %derivada de Ay em relação a x
            dAx_dy=(A(1,i,j+1,k)-A(1,i,j-1,k))/2/dy;%derivada de Ax em relação a y
            
            B(1,i-1,j-1,k-1) = dAz_dy-dAy_dz; %a componente x de J é dAz/dy-dAy/dz
            B(2,i-1,j-1,k-1) = dAx_dz-dAz_dx;%a componente y de J é dAx/dz-dAz/dx
            B(3,i-1,j-1,k-1) = dAy_dx-dAx_dy;%a componente z de J é dAy/dx-dAx/dy
        end
    end
end

xmedio3=ceil(length(xn2)/2);
ymedio3=ceil(length(yn2)/2);
zmedio3=ceil(length(zn2)/2);


%% Graficos B
%Gráfico XY de B em z = 0
figure(8);
[X,Y] = meshgrid(xn2,yn2);
quiver(X, Y, squeeze(B(1,:,:,zmedio3))', squeeze(B(2,:,:,zmedio3))');
hold on
contour(X, Y, squeeze(B(3,:,:,zmedio3))', 20);
% axis([-1 1 -1 1]);
xlabel('eixo x (m)')
ylabel('eixo y (m)')
title('B, plano XY (z=0)');

% 
%Gráfico XZ de M em y = 0
figure(6);
[X,Z] = meshgrid(xn2,zn2);
quiver(X, Z, squeeze(B(1,:,ymedio3,:))', squeeze(B(3,:,ymedio3,:))');
hold on
contour(X, Z, squeeze(B(2,:,ymedio3,:))', 20);
% axis([-1 1 -1 1]);
xlabel('eixo x (m)')
ylabel('eixo z (m)')
title('B, plano XZ (y=0)');