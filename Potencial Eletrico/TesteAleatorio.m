%% Descobrir o E1 - Campo eletrico dentro do cilindro
%tem que ser por gauss pq é um fio infinito
 

%% Questão 6
% No espaço livre, uma distribuição volumétrica de cargas constante ?v = 1 /m3
% existe dentro da região 0 ? ? ? 1 m, 0 ? ? ? 2? rad, e -1 m ? z ? 1 m. Calcule E,
% numericamente, em qualquer posição.Faça uma representação gráfica de E no
% plano xz.

%% Variaveis Dadas
clc
clear all
close all

%variaveis do problema
e0 = 8.854*10^-12;
pv = 8.1*10^-12;
a = 2.7*10^-3;
b = 7.3*10^-3;

k=1/(4*pi*e0);  %Constante

%% Variaveis Criadas

passe = 1/250;
limites = 1/10;
%Onde o campo sera medido:
x=[-limites:passe:limites]; %vetor na coordenada x onde será calculado E
y=[-limites:passe:limites]; %vetor na coordenada z onde será calculado E
%z=[-limites:passe:limites]; %vetor na coordenada z onde será calculado E

%Gerador do campo:
xl=[-passe*10:passe:passe*10]; % variação da coordenada x onde está a carga 
yl=[-passe*10:passe:passe*10];
zl=[-limites:passe:limites];

dV = passe; %área de cada segmento

%inicializa o campo elétrico:
V(1,:,:) = zeros (length(x),length(y)); 
V(2,:,:) = zeros (length(x),length(y));
V(3,:,:) = zeros (length(x),length(y));

%% Desenvolvimento
for i = 1:length(x)% varre a coordenada x onde E será calculado
    disp(i)
    for j = 1: length(y)  % varre a coordenada y onde E será calculado
        %for k = 1: length(z)  % varre a coordenada y onde E será calculado

            for m = 1:length(xl)  % varre a coordenada x da carga
                for n = 1:length(yl)  % varre a coordenada y da carga
                    for o = 1:length(zl)

                        r = [x(i),y(j),0]; %vetor posição apontando para onde estamos calculando E
                        rl= [xl(m),yl(n),zl(o)];% vetor posição apontando para um segmento do disco

                        if ((r-rl)*(r-rl)'>passe/10)
                            if (sqrt(xl(m)^2+yl(n)^2)<=a)
                                V(:,i,j) = V(:,i,j)  + k*((pv*dV)/sqrt((r-rl)*(r-rl)')); % para cada ponto (xl, yl)  do disco somo a contribuição da carga para o campo na posição (x,z).
                            
                            elseif (sqrt(xl(m)^2+yl(n)^2)>b)
                                V(:,i,j) = V(:,i,j)  + 0; % para cada ponto (xl, yl)  do disco somo a contribuição da carga para o campo na posição (x,z).
                            end
                        end

                    end
                end
            end
       % end
    end
end
%% Grafico
[X,Y] = meshgrid(x,y);

quiver(X,Y,squeeze(V(1,:,:))' , squeeze(V(2,:,:))')  %faz o gráfico vetorial
xlabel('eixo x (m)')
ylabel('eixo y (m)')


max = length(x)-1;
Vinf = V(1,max,max/2);
V0 =  V(1,max/2,max/2);
ddp = V0 - Vinf;
disp(ddp)
%0,62