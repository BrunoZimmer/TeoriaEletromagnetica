%% Descobrir o E1 - Campo eletrico dentro do cilindro


%% Variaveis Dadas
clc
clear all
close all

%% variaveis do problema
e0 = 8.854*10^-12;
kc = 1*10^-12;

const=1/(4*pi*e0);  %Constante

%% Variaveis Criadas

limites=30;
a=10;
dL=a/10; %tamanho dos segmentos

%Onde o campo sera medido:
x= [-limites:3:limites]; %varia��o da coordenada x onde ser� calculado V
y= [-limites:3:limites]; %varia��o da coordenada y onde ser� calculado V
z= [-limites:3:limites]; %varia��o da coordenada z onde ser� calculado V

midx = int64(length(x)/2);
midy = int64(length(y)/2);
midz = int64(length(z)/2);

%divis�o da linha em segmentos
xl=[a:dL:30*limites]; % varia��o da coordenada x onde est� a carga

V(:,:,:) = zeros (length(x), length(y), length(z));

%% Desenvolvimento
for i = 1: length(x)  % varre a coordenada x onde V ser� calculado
    disp(i)
    for j = 1: length(y)  % varre a coordenada y onde V ser� calculado
        
        for k = 1: length(z) % varre a coordenada z onde V ser� calculado
            
            for m = 1:length(xl)  % varre a coordenada x da carga
                
                r = [x(i),y(j),z(k)]; %vetor posi��o apontando para onde estamos calculando V
                rl= [xl(m),0,0];% vetor posi��o apontando para um segmento da linha
                
                %if ( (r-rl)*(r-rl)'>1/100)
                V(i,j,k) = V(i,j,k)  + const*((((kc.*xl(m))/(xl(m).^2 + a^2)).*dL)/(sqrt((r-rl)*(r-rl)'))');
                    
                %end
                
                
            end
        end
    end
end
%% Grafico
% xd = linspace(-limites,limites);
% yd = linspace(-limites,limites);
% zd = linspace(-limites,limites);
% [X,Y] = meshgrid(xd,yd);

[Y,Z] = meshgrid(y,z);
figure
[C,h] =  contour3(y,z,squeeze(V(midx,:,:)),20); %faz o gr�fico das curvas de n�vel para o potencial
set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
xlabel('eixo y (m)')
ylabel('eixo z (m)')
hold

[X,Z] = meshgrid(x,z);
figure
[C,h] = contour3(x,z,squeeze(V(:,midy,:)),20); %faz o gr�fico das curvas de n�vel para o potencial
set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
xlabel('eixo z (m)')
ylabel('eixo x (m)')


%% Print resultado
% max = int64(length(x));
% V0 =  V(max/2,max/2);
% Vinf =  0;
% disp(V0)