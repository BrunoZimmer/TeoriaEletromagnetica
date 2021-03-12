%% Quest�o 4
% Determine a intensidade de campo el�trico sobre o eixo z produzido por um anel
% de densidade superficial uniforme de carga ?s no espa�o livre. O anel ocupa a
% regi�o z=0, 1 m ? ? ? 2 m, 0 ? f ? 2? rad, em coordenadas cil�ndricas. Calcule,
% numericamente, a da intensidade de campo el�trico em qualquer posi��o,
% considerando ?s =1 pC/m2
% . Fa�a uma representa��o gr�fica de E no plano xz.

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
x=[-4:1/10:4]; %vetor na coordenada x onde ser� calculado E
z=[-4:1/10:4]; %vetor na coordenada z onde ser� calculado E

%Gerador do campo:
xl=[-2:1/10:2]; % varia��o da coordenada x onde est� a carga 
yl=[-2:1/10:2];

A = 0.01; %�rea de cada dA, area medida por vez

%inicializa o campo el�trico:
E(1,:,:) = zeros (length(x),length(z)); 
E(2,:,:) = zeros (length(x),length(z));
E(3,:,:) = zeros (length(x),length(z));

for i = 1:length(x)% varre a coordenada x onde E ser� calculado
    disp(i)%contador
    for j = 1: length(z)  % varre a coordenada y onde E ser� calculado
        for m = 1:length(xl)  % varre a coordenada x da carga
            for n = 1:length(yl)  % varre a coordenada y da carga
                
                r = [x(i),0,z(j)]; %vetor posi��o apontando para onde estamos calculando E
                rl= [xl(m),yl(n),0];% vetor posi��o apontando para um segmento do disco
                
                %condi��o para se ter os dois "circulos" do anel(r=1 e r=2)
                %e n�o dar erro numero de proximidade
                if (sqrt(xl(m)^2+yl(n)^2)<=2) && (sqrt(xl(m)^2+yl(n)^2)>=1) && ((r-rl)*(r-rl)'>1/100)
                    E(:,i,j) = E(:,i,j)  + (ps*A/sqrt((r-rl)*(r-rl)')^3*(r-rl))'; % para cada ponto (xl, yl)  do disco somo a contribui��o da carga para o campo na posi��o (x,z).
                end
            end
        end
    end
end
E= E*k;%formula de coulomb

[X,Z] = meshgrid(x,z);
quiver(X,Z,squeeze(E(1,:,:))' , squeeze(E(3,:,:))')  %faz o gr�fico vetorial
xlabel('eixo x (m)')
ylabel('eixo z (m)')

