%Q3
% Considere uma linha de carga, no eixo z, com densidade linear de cargas
% uniforme pL = 1 C/m, se estendendo de z = -1 m at� z = 1 m. Calcule,
% numericamente, a intensidade de campo el�trico E em qualquer posi��o espa�o.
% Fa�a uma representa��o gr�fica de E no plano xz e uma no plano xy.


clc
clear all
close all

z1=-1;
z2=1;
pL=1;  %Densidade Superficial de Carga [C/m^2]
k=1/(4*pi*8.854*10^(-12));  %Constante

x=[-2:1/20:2]; %vetor na coordenada x onde ser� calculado E
z=[-2:1/20:2]; %vetor na coordenada z onde ser� calculado E

zl=[z1:1/20:z2]; %  varia��o da coordenada x onde est� a carga 
A = 1; %�rea de cada dA, area medida por vez

%inicializa o campo el�trico:
E(1,:,:) = zeros (length(x),length(z));
E(2,:,:) = zeros (length(x),length(z));
E(3,:,:) = zeros (length(x),length(z));

for i = 1:length(x)% varre a coordenada x onde E ser� calculado
    disp(i)%contador
    for j = 1: length(z)  % varre a coordenada z onde E ser� calculado
        for m = 1:length(zl)  % varre a coordenada x da carga
            
            r = [x(i),0,z(j)]; %vetor posi��o apontando para onde estamos calculando E
            rl= [0,0,zl(m)];% vetor posi��o apontando para um segmento da placa

            if ((r-rl)*(r-rl)'>1/100)%nao dar erro numerico de coisas mt proximas
                E(:,i,j) = E(:,i,j)  + (pL*A/sqrt((r-rl)*(r-rl)')^3*(r-rl))'; % para cada ponto (xl, yl)  da placa somo a contribui��o da carga para o campo na posi��o (x,z). Considero a carga concentrada no centro do segmento.
            end
        end
    end
end
E= E*k;%formula de coulomb

[X,Z] = meshgrid(x,z);
quiver(X,Z,squeeze(E(1,:,:))' , squeeze(E(3,:,:))')  %faz o gr�fico vetorial
xlabel('eixo x (m)')
ylabel('eixo z (m)')

