%% Questão 4
% Determine a intensidade de campo elétrico sobre o eixo z produzido por um anel
% de densidade superficial uniforme de carga ?s no espaço livre. O anel ocupa a
% região z=0, 1 m ? ? ? 2 m, 0 ? f ? 2? rad, em coordenadas cilíndricas. Calcule,
% numericamente, a da intensidade de campo elétrico em qualquer posição,
% considerando ?s =1 pC/m2
% . Faça uma representação gráfica de E no plano xz.

%% Variaveis Dadas
clc
clear all
close all

x1=-2;
x2=-1;
x3=1;
x4=2;

ps=1*10^(-12);  %Densidade Superficial de Carga [C/m^2]
k=1/(4*pi*8.854*10^(-12));  %Constante


%Onde o campo sera medido:
x=[-4:1/10:4]; %vetor na coordenada x onde será calculado E
z=[-4:1/10:4]; %vetor na coordenada z onde será calculado E

%Gerador do campo:
xl=[-2:1/10:2]; % variação da coordenada x onde está a carga 
yl=[-2:1/10:2];

A = 0.01; %área de cada dA, area medida por vez

%inicializa o campo elétrico:
E(1,:,:) = zeros (length(x),length(z)); 
E(2,:,:) = zeros (length(x),length(z));
E(3,:,:) = zeros (length(x),length(z));

for i = 1:length(x)% varre a coordenada x onde E será calculado
    disp(i)%contador
    for j = 1: length(z)  % varre a coordenada y onde E será calculado
        for m = 1:length(xl)  % varre a coordenada x da carga
            for n = 1:length(yl)  % varre a coordenada y da carga
                
                r = [x(i),0,z(j)]; %vetor posição apontando para onde estamos calculando E
                rl= [xl(m),yl(n),0];% vetor posição apontando para um segmento do disco
                
                %condição para se ter os dois "circulos" do anel(r=1 e r=2)
                %e não dar erro numero de proximidade
                if (sqrt(xl(m)^2+yl(n)^2)<=2) && (sqrt(xl(m)^2+yl(n)^2)>=1) && ((r-rl)*(r-rl)'>1/100)
                    E(:,i,j) = E(:,i,j)  + (ps*A/sqrt((r-rl)*(r-rl)')^3*(r-rl))'; % para cada ponto (xl, yl)  do disco somo a contribuição da carga para o campo na posição (x,z).
                end
            end
        end
    end
end
E= E*k;%formula de coulomb

[X,Z] = meshgrid(x,z);
quiver(X,Z,squeeze(E(1,:,:))' , squeeze(E(3,:,:))')  %faz o gráfico vetorial
xlabel('eixo x (m)')
ylabel('eixo z (m)')

