% Uma espira filamentar quadrada de lado a está inicialmente posicionada
% no plano xy, com seu centro alinhado com o eixo z, numa região de campo 
% uniforme B(y,t) = âz B0 cos(?t) + ây B0 sen(2?t) . A espira começa a 
% girar com velocidade angular ? sobre o eixo y. Determine numericamente
% o fluxo magnético através da espira e a FEM induzida nessa espira ao 
% longo do tempo. Trace os gráficos das duas quantidades em função do tempo.




clear all
clc
close all

%% Variaveis dadas

u0=4*pi*10^(-7); % Permeabilidade magnÃ©tica do espaÃ§o livre (em H/m)
a = 1;
B0 = 8;
w = 55;

tmp = 1;


%% Espaços de calculos

passo=a/5; % passo

dx=passo;
dy=passo;
dz=passo;

x= -2*a:dx:2*a; %variacao da coordenada x
y= -2*a:dy:2*a; %variacao da coordenada y
z= -2*a:dz:2*a; %variacao da coordenada z

xl= -a/2:dx:a/2;
yl= -a/2:dz:a/2;

xmedio = ceil(length(x)/2);
ymedio = ceil(length(y)/2);
zmedio = ceil(length(z)/2);


f=w/(2*pi);
T=1/f;

dt = T/100;
t = 0:dt:2*tmp;

%% Descobrindo valor do campo B

%zerando a densidade de fluxo magnético
B(:,:,:,:,:) = zeros(3, length(x), length(y), length(z), length(t));

for p = 1:length(t)
    for i = 1:length(x)
        disp(i)
        for j = 1:length(y)  
            for k = 1:length(z) 
                %B(y,t) = âz B0 cos(?t) + ây B0sen(2?t) 
                B(2,i,j,k,p) = B0*sin( 2*w*t(p) );
                B(3,i,j,k,p) = B0*cos( w*t(p) );
                
            end
        end
    end
    
end



%% Descobrindo o fluxo

fluxo = zeros(1, length(t));
for i = 1:length(xl)
    disp(i)
    for j = 1:length(yl)
        for p = 1:length(t)
            
            dS=dx*dy;
            n = [-sin(w*t(p)),0,cos(w*t(p))];
            
            dSn=n*dS;
            
            fluxo(1,p) = fluxo(1,p) + dot(B(:,i,j,zmedio),dSn);
        end
    end
end


%% Plot Flx x tmp
figure(1);
plot(t,fluxo);
xlabel('Tempo')
ylabel('W')
title('Flx x tmp');

%% Descobrindo a FEM
FEM = zeros(1,length(t));

for i = 2:length(t)-1% varre a coordenada x onde E será calculado
    disp(i);
    FEM(1,i) = (-1 * (fluxo(1,i+1)-fluxo(1,i-1))) / (2*dt);
end

%% Plot FEM x tmp
figure(2);
plot(t,FEM);
xlabel('Tempo')
ylabel('V')
title('FEM x tmp');

%% Disp dos valores de teste
disp(fluxo(ceil(0.5*length(t))));
disp(FEM(ceil(0.5*length(t))));


%% Plot animado do campo B em YZ

for i = 2:length(t)
    figure(3)
    [Y,Z] = meshgrid(y,z);
    h = quiver(Y, Z, squeeze(B(2,xmedio,:,:,i))', squeeze(B(3,xmedio,:,:,i))');
    xlabel('Eixo y')
    ylabel('Eixo z')
    drawnow 
    pause(0.1) 
end

for i = 2:length(t)
    figure(3)
    [X,Y,Z] = meshgrid(x,y,z);
    quiver3( X, Y, Z, squeeze(B(1,:,:,:,i)),squeeze(B(2,:,:,:,i)),squeeze(B(3,:,:,:,i)));
%     AXIS([XMIN XMAX YMIN YMAX ZMIN ZMAX])
    axis([-0.8 0.8 -0.8 0.8 -0.8 0.8] )
    drawnow 
    pause(0.1) 
end