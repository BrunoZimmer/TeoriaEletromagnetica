%% Questao 1
% Considere uma placa retangular de lados L =1 m e W = 1 m, com espessura 
% desprezível.
% Nessa placa as cargas elétricas estão distribuídas uniformemente, com 
% ps = 1 C/m2.
% A placa está posicionada em z = 0 e centrada no eixo z. Calcule, 
% numericamente, a intensidade de campo elétrico E em qualquer posição do 
% plano xz. Faça uma representação gráfica de E.

%% Variaveis Dadas
clc
clear all
close all

L=1;    %Dimensão [metros]
W=1;    %Dimensão [metros]
ps=1;  %Densidade Superficial de Carga [C/m^2]
k=1/(4*pi*8.854*10^(-12));  %Constante

%% Variaveis Criadas
D=max(L,W);
%Onde o campo sera medido:
x=[-2*D:D/20:2*D]; %vetor na coordenada x onde será calculado E
z=[-3*D:D/20:3*D]; %vetor na coordenada z onde será calculado E

%Gerador do campo:
%vou dividir a placa em segmentos de dimensão L/20 x W/20
xl=[-L/2+L/40:L/20:L/2-L/40]; % variação da coordenada x onde está a carga (centro dos segmentos)
yl=[-W/2+W/40:W/20:W/2-W/40]; % variação da coordenada y onde está a carga
A = L/20 * W/20; %área de cada segmento

%inicializa o campo elétrico:
E(1,:,:) = zeros (length(x),length(z)); 
E(2,:,:) = zeros (length(x),length(z));
E(3,:,:) = zeros (length(x),length(z));

%% Desenvolvimento
for i = 1:length(x)% varre a coordenada x onde E será calculado
    i
    for j = 1: length(z)  % varre a coordenada z onde E será calculado
        
        for m = 1:length(xl)  % varre a coordenada x da carga
            for n=1:length(yl) % varre a coordenada  y da carga
                
                r = [x(i),0,z(j)]; %vetor posição apontando para onde estamos calculando E
                rl= [xl(m),yl(n),0];% vetor posição apontando para um segmento da placa
                
                %if ( (r-rl)*(r-rl)'>D/100)
                E(:,i,j) = E(:,i,j)  + (ps*A/sqrt((r-rl)*(r-rl)')^3*(r-rl))'; % para cada ponto (xl, yl)  da placa somo a contribuição da carga para o campo na posição (x,z). Considero a carga concentrada no centro do segmento.
                %end

            end
        end
        
    end
end
E= E/k;

%% Grafico
[X,Z] = meshgrid(x,z);
quiver(X,Z,squeeze(E(1,:,:))' , squeeze(E(3,:,:))')  %faz o gráfico vetorial
xlabel('eixo x (m)')
ylabel('eixo z (m)')

