%-------------------------------------------------------------------------%
%A superf�cie anelar 1 cm < ? < 3 cm, z =0, est� carregada com a densidade 
%superficial n�o uniforme de cargas ps = 6,8? nC/m2. Assumindo V=0 no infinito,
%determine numericamente o potencial devido � superficie anelar carregada em 
%qualquer ponto do espa�o. Com os dados obtidos, trace o gr�fico do potencial
%em fun��o de z, considerando x=0 e y=0. Trace tamb�m o gr�fico do potencial 
%nos planos xz e xy. Para validar sua resposta, calcule V em P(x = 0, y = 0, z= 2) cm.
%-------------------------------------------------------------------------%

clc
close all; 
clear all;
%% Constantes do problema
ri=0.01;
rf=0.03;

const = 1/(4*pi*8.854*10^(-12));
ps = 7.9*10^-12; 

%% Limites
limites=1.4*rf;
passe=ri/4; %tamanho dos segmentos 

%Onde o campo sera medido:
x= -limites:passe:limites; %vetor na coordenada x onde ser� calculado V
y= -limites:passe:limites;  %vetor na coordenada z onde ser� calculado V
z= -limites:passe:limites; %vetor na coordenada z onde ser� calculado V

%Gerador do campo:
xl=[-rf:passe/3:rf]; % varia��o da coordenada x onde est� a carga
yl=[-rf:passe/3:rf]; % varia��o da coordenada y onde est� a carga

dS=passe^2; %tamanho de cada segmento

%inicializa o campo:
V(:,:,:) = zeros (length(x), length(y), length(z));

%% Desenvolvimento
for i = 1: length(x)  % varre a coordenada y onde V ser� calculado
    disp(i);
    for j = 1: length(y)  % varre a coordenada y onde V ser� calculado
        for k = 1: length(z) % varre a coordenada z onde V ser� calculado
            
            for m = 1:length(xl)  % varre a coordenada x da carga
                for n = 1:length(yl)  % varre a coordenada y da carga
                    
                    r = [x(i),y(j),z(k)]; %vetor posi��o apontando para onde estamos calculando E
                    rl= [xl(m),yl(n),0];% vetor posi��o apontando para um segmento da linha
                    
%                         if ((r-rl)*(r-rl)'>1/100)
                    if (sqrt(xl(m)^2+yl(n)^2) > ri) && (sqrt(xl(m)^2+yl(n)^2) <= rf) 
                        V(i,j,k) = V(i,j,k)  + (const*(ps*sqrt(xl(m)^2+yl(n)^2)*dS)/sqrt((r-rl)*(r-rl)'))';
                    end
%                         end       
                     
                 end
            end
            
        end
    end
end

mid = int64(length(x)/2);

[X,Y] = meshgrid(x,y);
figure
[C,h] =  contour(x,y,squeeze(V(:,:,mid)),30); %faz o gr�fico das curvas de n�vel para o potencial
%set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
xlabel('eixo x (m)')
ylabel('eixo y (m)')
hold

[X,Y] = meshgrid(x,y);
figure
contour3(x,y,squeeze(V(:,:,mid)),30); %faz o gr�fico das curvas de n�vel para o potencial
%set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
xlabel('eixo x (m)')
ylabel('eixo y (m)')
hold

[X,Z] = meshgrid(x,z);
figure
[C,h] =  contour(x,z,squeeze(V(:,mid,:)),30); %faz o gr�fico das curvas de n�vel para o potencial
%set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
xlabel('eixo x (m)')
ylabel('eixo z (m)')
hold
 