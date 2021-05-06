% Por um pedaço de fio dobrado na forma de espira, conforme a figura abaixo,
% passa uma corrente que aumenta com o tempo I(t) = kt (-? < t < ?) Calcule
% o potencial vetorial retardado A em qualquer posição do espaço. Calcule B
% a partir desse resultado. Encontre a intensidade de campo elétrico E  em
% qualquer posição.
% 
% Faça representações gráficas de A, E e B nos planos xy, xz e yz. 
% 
% Afim de avaliar a sua resposta calcule a magnitude da intensidade de
% campo elétrico , na posição x, y e z =0, considerando k = 8,2 A/s,
% b = 5,0 m e a = 2,3 m

%% Variaveis iniciais

u0=4*pi*10^(-7); 
c = 3*10^(8);
k = 8.2;
b = 5.2;
a = 2;
tmp = 0.5;

%% Espaços de calculos

passo=a/3; % passo

dx=passo;
dy=passo;
dz=passo;

x= -1.5*b:dx:1.5*b; %variacao da coordenada x
y= -1.5*b:dy:1.5*b; %variacao da coordenada y
z= -1.5*b:dz:1.5*b; %variacao da coordenada z

xmedio = ceil(length(x)/2);
ymedio = ceil(length(y)/2);
zmedio = ceil(length(z)/2);

dphi = pi/10;
%parte 1 da espira / semi-circulo a 
pl1= a;
% phil1= 0:dphi:pi;
phil1= 0:dphi:pi;

%parte 2 da espira / semi-circulo b
pl2= b;
phil2= 0:dphi:pi;

%parte 3 da espira / reta de a a b
xl1 = a:dx:b;
yl1 = 0;

%parte 4 da espira / reta de -b a -a
xl2 = -b:dx:-a;
yl2 = 0;

%posicao da espira no eixo z
h=0;

%variação de tmepo
dt = tmp/20;
t = 0:dt:2*tmp;
tmedio  = ceil(length(t)/2);

%% Descobrindo valor do campo A

%zerando a densidade de fluxo magnético
A = zeros(3, length(x), length(y), length(z), length(t));

for p = 1:length(t)
    for i = 1:length(x)
        disp(i)
        for j = 1:length(y)  
            for g = 1:length(z) 
                
                r = [x(i),y(j),z(g)];  %vetor posiÃ§Ã£o do ponto de campo (onde calculamos B)
                
%                 A = int(I/r) = int(kt/r)  =  int(k(t-r/c)/r)
                for u = 1:length(phil1)  
                    [rl1(1), rl1(2), rl1(3)] = pol2cart(phil1(u), pl1, h); % convertemos rl para coord retangulares

                    [Ir1(1),Ir1(2),Ir1(3)] = pol2cart(phil1(u)-pi/2, (k*(t(p)- norm(r-rl1)/c)) , h);
                    dA1 = (Ir1 /norm(r-rl1));
                    A(1,i,j,g,p) = A(1,i,j,g,p) + dA1(1)*a*dphi;
                    A(2,i,j,g,p) = A(2,i,j,g,p) + dA1(2)*a*dphi;
                    A(3,i,j,g,p) = A(3,i,j,g,p) + dA1(3)*a*dphi;
                end
                
%                 norm(r-rl2) = raiz(r² - rl²)-> a distancia da carga até o
%                 ponto de analise

% 
             
                
                
                for v = 1:length(phil2)  
                    [rl2(1), rl2(2), rl2(3)] = pol2cart(phil2(v), pl2, h); % convertemos rl para coord retangulares
                    
                    [Ir2(1),Ir2(2),Ir2(3)] = pol2cart(phil2(v)+pi/2, (k*(t(p)- norm(r-rl2)/c)) , h);
%                     if  ((r-rl2)*(r-rl2)' > (dphi/2))
                    dA2 = (Ir2 /norm(r-rl2));
                    A(1,i,j,g,p) = A(1,i,j,g,p) + dA2(1)*a*dphi;
                    A(2,i,j,g,p) = A(2,i,j,g,p) + dA2(2)*a*dphi;
                    A(3,i,j,g,p) = A(3,i,j,g,p) + dA2(3)*a*dphi;
                end
                
                
                for w = 1:length(xl1)  
                    rl3 =[xl1(w), yl1,h]; % convertemos rl para coord retangulares
                    
                    Ir3 = [(k*(t(p)- norm(r-rl3)/c)), 0, 0];
                    
%                     if  ((r-rl3)*(r-rl3)' > (dphi/2))
                    dA3 = (Ir3 /norm(r-rl3));

                    A(1,i,j,g,p) = A(1,i,j,g,p) + dA3(1)*dx;
                    A(2,i,j,g,p) = A(2,i,j,g,p) + dA3(2)*dx;
                    A(3,i,j,g,p) = A(3,i,j,g,p) + dA3(3)*dx;
                end
                
                for f = 1:length(xl2)  
                    rl4 =[xl2(f), yl2,h]; % convertemos rl para coord retangulares
                    
                    Ir4 = [(k*(t(p)- norm(r-rl4)/c)), 0, 0];
                    
                    dA4 = (Ir4 /norm(r-rl4));

                    A(1,i,j,g,p) = A(1,i,j,g,p) + dA4(1)*dx;
                    A(2,i,j,g,p) = A(2,i,j,g,p) + dA4(2)*dx;
                    A(3,i,j,g,p) = A(3,i,j,g,p) + dA4(3)*dx;
                end
            end
        end
    end
end

A = ( u0/(4*pi) ) * A;

%% Plot A

%gráfico vetorial de A no plano XY
[X,Y] = meshgrid(x,y);

figure (1);
quiver(X, Y, squeeze(A(1,:,:,zmedio,tmedio))', squeeze(A(2,:,:,zmedio,tmedio))'); 
% axis([-0.004 0.004 0.022 0.032])
xlabel('eixo x (m)')
ylabel('eixo y (m)')
title('A, plano XY (z=0)');

%gráfico vetorial de A no plano YZ
[Y,Z] = meshgrid(y,z);

figure (2);
quiver(Y, Z, squeeze(A(2,xmedio,:,:,tmedio))', squeeze(A(3,xmedio,:,:,tmedio))'); 
xlabel('eixo y (m)')
ylabel('eixo z (m)')
title('A, plano YZ (x=0)');

%gráfico vetorial de A no plano YZ
[X,Z] = meshgrid(x,z);

figure (3);
quiver(X, Z, squeeze(A(2,:,ymedio,:,tmedio))', squeeze(A(3,:,ymedio,:,tmedio))'); 
xlabel('eixo x (m)')
ylabel('eixo z (m)')
title('A, plano XZ (y=0)');

%% Calculo de B a partir do rotacional de A

% rotacional de A
% B = DELTA x A

%Derivada é um espaço menor tirando o "cantos" do espaço
xn = (x(2:end-1));
yn = (y(2:end-1));
zn = (z(2:end-1));
tn = (t(2:end-1));

%Zerando densidades de corrente ligadas
B  = zeros (3,length(xn),length(yn),length(zn),length(t)); 

%rot M = (dAz/dy-dAy/dz) ax + (dAx/dz-dAz/dx) ay + (dAy/dx-dAx/dy) az
for p = 1:length(t)
    for i = 2:length(x)-1
        disp(i)
        for j = 2 :length(y)-1 
            for k = 2: length(z)-1  

                dAz_dy= (A(3,i,j+1,k,p)-A(3,i,j-1,k,p))/2/dy; %derivada de Az em relação a y
                dAy_dz = (A(2,i,j,k+1,p)-A(2,i,j,k-1,p))/2/dz; %derivada de Ay em relação a z

                dAx_dz=(A(1,i,j,k+1,p)-A(1,i,j,k-1,p))/2/dz; %derivada de Ax em relação a z
                dAz_dx=(A(3,i+1,j,k,p)-A(3,i-1,j,k,p))/2/dx;%derivada de Az em relação a x

                dAy_dx=(A(2,i+1,j,k,p)-A(2,i-1,j,k,p))/2/dx; %derivada de Ay em relação a x
                dAx_dy=(A(1,i,j+1,k,p)-A(1,i,j-1,k,p))/2/dy;%derivada de Ax em relação a y

                B(1,i-1,j-1,k-1,p) = dAz_dy-dAy_dz; %a componente x de B é dAz_dy-dAy_dz
                B(2,i-1,j-1,k-1,p) = dAx_dz-dAz_dx;%a componente y de B é dAx_dz-dAz_dx
                B(3,i-1,j-1,k-1,p) = dAy_dx-dAx_dy;%a componente z de B é dAy_dx-dAx_dy

            end
        end
    end
end


%% Plot B

xmedio2=ceil(length(xn)/2);
ymedio2=ceil(length(yn)/2);
zmedio2=ceil(length(zn)/2);
tmedio2=ceil(length(tn)/2);

%gráfico vetorial de B no plano XY
[X,Y] = meshgrid(xn,yn);

figure (4);
quiver(X, Y, squeeze(B(1,:,:,zmedio2,tmedio))', squeeze(B(2,:,:,zmedio2,tmedio))');
hold on
contour(xn, yn, squeeze(B(3,:,:,zmedio2,tmedio))', 50); 
% axis([-0.004 0.004 0.022 0.032])
xlabel('eixo x (m)')
ylabel('eixo y (m)')
title('B, plano XY (z=0)');

%gráfico vetorial de B no plano YZ
[Y,Z] = meshgrid(yn,zn);

figure (5);
quiver(Y, Z, squeeze(B(2,xmedio2,:,:,tmedio))', squeeze(B(3,xmedio2,:,:,tmedio))'); 
xlabel('eixo y (m)')
ylabel('eixo z (m)')
title('B, plano YZ (x=0)');

%gráfico vetorial de B no plano YZ
[X,Z] = meshgrid(xn,zn);

figure (6);
quiver(X, Z, squeeze(B(2,:,ymedio2,:,tmedio))', squeeze(B(3,:,ymedio2,:,tmedio))'); 
xlabel('eixo x (m)')
ylabel('eixo z (m)')
title('B, plano XZ (y=0)');


%% Calculo de E a partir da derivada negativa de A

%Zerando densidades de corrente ligadas
E = zeros (3,length(x),length(y),length(z),length(tn)); 

%E= -dA/dt
for p = 2:length(t)-1
    for i = 2:length(x)-1
        disp(i)
        for j = 2 :length(y)-1 
            for k = 2: length(z)-1  

                dAx_dt= (A(1,i,j,k,p+1)-A(1,i,j,k,p-1))/(2*dt); %derivada de Az em relação a t
                dAy_dt= (A(2,i,j,k,p+1)-A(2,i,j,k,p-1))/(2*dt); %derivada de Az em relação a t
                dAz_dt= (A(3,i,j,k,p+1)-A(3,i,j,k,p-1))/(2*dt); %derivada de Az em relação a t

                E(1,i,j,k,p-1) = dAx_dt;
                E(2,i,j,k,p-1) = dAy_dt;
                E(3,i,j,k,p-1) = dAz_dt;

            end
        end
    end
end


%% Plot E

%gráfico vetorial de E no plano XY
[X,Y] = meshgrid(x,y);

figure (7);
quiver(X, Y, squeeze(E(1,:,:,zmedio,tmedio2))', squeeze(E(2,:,:,zmedio2,tmedio2))'); 
% axis([-0.004 0.004 0.022 0.032])
xlabel('eixo x (m)')
ylabel('eixo y (m)')
title('E, plano XY (z=0)');

%gráfico vetorial de E no plano YZ
[Y,Z] = meshgrid(y,z);

figure (8);
quiver(Y, Z, squeeze(E(2,xmedio,:,:,tmedio2))', squeeze(E(3,xmedio,:,:,tmedio2))'); 
xlabel('eixo y (m)')
ylabel('eixo z (m)')
title('E, plano YZ (x=0)');

%gráfico vetorial de E no plano YZ
[X,Z] = meshgrid(x,z);

figure (9);
quiver(X, Z, squeeze(E(2,:,ymedio,:,tmedio2))', squeeze(E(3,:,ymedio,:,tmedio2))'); 
xlabel('eixo x (m)')
ylabel('eixo z (m)')
title('E, plano XZ (y=0)');

