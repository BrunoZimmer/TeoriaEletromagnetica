%% Questão 6
% No espaço livre, uma distribuição volumétrica de cargas constante ?v = 1 /m3
% existe dentro da região 0 ? ? ? 1 m, 0 ? ? ? 2? rad, e -1 m ? z ? 1 m. Calcule E,
% numericamente, em qualquer posição.Faça uma representação gráfica de E no
% plano xz.

%% Variaveis Dadas
clc
clear all
close all

x1=-1;
x2=1;
y1=-1;
y2=1;
z1=-1;
z2=1;

pv=1*10^(-12);  %Densidade Superficial de Carga [C/m^2]
k=1/(4*pi*8.854*10^(-12));  %Constante

%% Variaveis Criadas

%Onde o campo sera medido:
x=[-2:1/10:2]; %vetor na coordenada x onde será calculado E
z=[-2:1/10:2]; %vetor na coordenada z onde será calculado E

%Gerador do campo:
xl=[-1:1/10:1]; % variação da coordenada x onde está a carga 
yl=[-1:1/10:1];
zl=[-1:1/10:1];

V = 0.1; %área de cada segmento

%inicializa o campo elétrico:
E(1,:,:) = zeros (length(x),length(z)); 
E(2,:,:) = zeros (length(x),length(z));
E(3,:,:) = zeros (length(x),length(z));

%% Desenvolvimento
for i = 1:length(x)% varre a coordenada x onde E será calculado
    disp(i)
    for j = 1: length(z)  % varre a coordenada y onde E será calculado
        
        for m = 1:length(xl)  % varre a coordenada x da carga
            for n = 1:length(yl)  % varre a coordenada y da carga
                for o = 1:length(zl)
                    
                    r = [x(i),0,z(j)]; %vetor posição apontando para onde estamos calculando E
                    rl= [xl(m),yl(n),zl(o)];% vetor posição apontando para um segmento do disco

                    if (sqrt(xl(m)^2+yl(n)^2)<=1)
                        if ((r-rl)*(r-rl)'>1/100)
                            E(:,i,j) = E(:,i,j)  + (pv*V/sqrt((r-rl)*(r-rl)')^3*(r-rl))'; % para cada ponto (xl, yl)  do disco somo a contribuição da carga para o campo na posição (x,z).
                        end
                    end
                    
                end
            end
        end
    end
end
E= E*k;%formula de coulomb

%% Grafico
[X,Z] = meshgrid(x,z);
quiver(X,Z,squeeze(E(1,:,:))' , squeeze(E(3,:,:))')  %faz o gráfico vetorial
xlabel('eixo x (m)')
ylabel('eixo z (m)')
