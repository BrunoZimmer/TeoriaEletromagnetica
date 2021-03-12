%% Questao 1
% Considere uma placa retangular de lados L =1 m e W = 1 m, com espessura 
% desprez�vel.
% Nessa placa as cargas el�tricas est�o distribu�das uniformemente, com 
% ps = 1 C/m2.
% A placa est� posicionada em z = 0 e centrada no eixo z. Calcule, 
% numericamente, a intensidade de campo el�trico E em qualquer posi��o do 
% plano xz. Fa�a uma representa��o gr�fica de E.

%% Variaveis Dadas
clc
clear all
close all

L=1;    %Dimens�o [metros]
W=1;    %Dimens�o [metros]
ps=1;  %Densidade Superficial de Carga [C/m^2]
k=1/(4*pi*8.854*10^(-12));  %Constante

%% Variaveis Criadas
D=max(L,W);
%Onde o campo sera medido:
x=[-2*D:D/20:2*D]; %vetor na coordenada x onde ser� calculado E
z=[-3*D:D/20:3*D]; %vetor na coordenada z onde ser� calculado E

%Gerador do campo:
%vou dividir a placa em segmentos de dimens�o L/20 x W/20
xl=[-L/2+L/40:L/20:L/2-L/40]; % varia��o da coordenada x onde est� a carga (centro dos segmentos)
yl=[-W/2+W/40:W/20:W/2-W/40]; % varia��o da coordenada y onde est� a carga
A = L/20 * W/20; %�rea de cada segmento

%inicializa o campo el�trico:
E(1,:,:) = zeros (length(x),length(z)); 
E(2,:,:) = zeros (length(x),length(z));
E(3,:,:) = zeros (length(x),length(z));

%% Desenvolvimento
for i = 1:length(x)% varre a coordenada x onde E ser� calculado
    i
    for j = 1: length(z)  % varre a coordenada z onde E ser� calculado
        
        for m = 1:length(xl)  % varre a coordenada x da carga
            for n=1:length(yl) % varre a coordenada  y da carga
                
                r = [x(i),0,z(j)]; %vetor posi��o apontando para onde estamos calculando E
                rl= [xl(m),yl(n),0];% vetor posi��o apontando para um segmento da placa
                
                %if ( (r-rl)*(r-rl)'>D/100)
                E(:,i,j) = E(:,i,j)  + (ps*A/sqrt((r-rl)*(r-rl)')^3*(r-rl))'; % para cada ponto (xl, yl)  da placa somo a contribui��o da carga para o campo na posi��o (x,z). Considero a carga concentrada no centro do segmento.
                %end

            end
        end
        
    end
end
E= E/k;

%% Grafico
[X,Z] = meshgrid(x,z);
quiver(X,Z,squeeze(E(1,:,:))' , squeeze(E(3,:,:))')  %faz o gr�fico vetorial
xlabel('eixo x (m)')
ylabel('eixo z (m)')

