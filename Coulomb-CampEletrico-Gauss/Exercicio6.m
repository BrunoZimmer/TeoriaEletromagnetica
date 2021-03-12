%% Quest�o 6
% No espa�o livre, uma distribui��o volum�trica de cargas constante ?v = 1 /m3
% existe dentro da regi�o 0 ? ? ? 1 m, 0 ? ? ? 2? rad, e -1 m ? z ? 1 m. Calcule E,
% numericamente, em qualquer posi��o.Fa�a uma representa��o gr�fica de E no
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
x=[-2:1/10:2]; %vetor na coordenada x onde ser� calculado E
z=[-2:1/10:2]; %vetor na coordenada z onde ser� calculado E

%Gerador do campo:
xl=[-1:1/10:1]; % varia��o da coordenada x onde est� a carga 
yl=[-1:1/10:1];
zl=[-1:1/10:1];

V = 0.1; %�rea de cada segmento

%inicializa o campo el�trico:
E(1,:,:) = zeros (length(x),length(z)); 
E(2,:,:) = zeros (length(x),length(z));
E(3,:,:) = zeros (length(x),length(z));

%% Desenvolvimento
for i = 1:length(x)% varre a coordenada x onde E ser� calculado
    disp(i)
    for j = 1: length(z)  % varre a coordenada y onde E ser� calculado
        
        for m = 1:length(xl)  % varre a coordenada x da carga
            for n = 1:length(yl)  % varre a coordenada y da carga
                for o = 1:length(zl)
                    
                    r = [x(i),0,z(j)]; %vetor posi��o apontando para onde estamos calculando E
                    rl= [xl(m),yl(n),zl(o)];% vetor posi��o apontando para um segmento do disco

                    if (sqrt(xl(m)^2+yl(n)^2)<=1)
                        if ((r-rl)*(r-rl)'>1/100)
                            E(:,i,j) = E(:,i,j)  + (pv*V/sqrt((r-rl)*(r-rl)')^3*(r-rl))'; % para cada ponto (xl, yl)  do disco somo a contribui��o da carga para o campo na posi��o (x,z).
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
quiver(X,Z,squeeze(E(1,:,:))' , squeeze(E(3,:,:))')  %faz o gr�fico vetorial
xlabel('eixo x (m)')
ylabel('eixo z (m)')
