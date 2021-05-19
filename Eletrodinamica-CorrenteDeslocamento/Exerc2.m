% Uma espira filamentar quadrada de lado a está inicialmente posicionada
% no plano xy, com seu centro alinhado com o eixo z, numa região de campo
% uniforme B(y,t) = B0 az. A espira começa a girar com velocidade angular
% ? sobre o eixo x. Determine numericamente o fluxo magnético através da
% espira e a FEM induzida nessa espira ao longo do tempo. Trace os gráficos
% das duas quantidades em função do tempo.

clear all
clc
close all

%% Variaveis dadas

u0=4*pi*10^(-7); % Permeabilidade magnÃ©tica do espaÃ§o livre (em H/m)
% a = 1;
% B0 = 8.1;
% w = 54.9;
% beta = 2.4;
% tmp = 2.8;

a = 0.05;        %lado da espira
B0 = 4.5;        %campo uniforme
w = 57;          %velocidade angular que a espira gira em torno do eixo x

tmp = 9.9;

%% Espaços de calculos

passo=a/10; % passo

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

dt = T/200;
t = 0:dt:2*tmp;

%% Descobrindo valor do campo B

%zerando a densidade de fluxo magnético
B(:,:,:,:) = zeros(3, length(x), length(y), length(z));
B(3,:,:,:) = B0;

%% Descobrindo o fluxo

fluxo = zeros(1, length(t));
for i = 1:length(xl)
    disp(i)
    for j = 1:length(yl)
        for p = 1:length(t)
            
            dS=dx*dy;
            n = [0,-sin(w*t(p)),cos(w*t(p))];
            
            dSn=n*dS;
            
            fluxo(1,p) = fluxo(1,p) + dot(B(:,i,j,zmedio),dSn);
        end
    end
end


%% Plot Flx x tmp
figure(2);
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
figure(3);
plot(t,FEM);
xlabel('Tempo')
ylabel('V')
title('FEM x tmp');

%% Disp dos valores de teste
disp(fluxo(ceil(0.5*length(t))));
disp(FEM(ceil(0.5*length(t))));





