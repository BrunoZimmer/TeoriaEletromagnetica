%% Questão 7
% Na região do espaço livre que inclui o volume 1,8 m < x < 5 m, 1,8 m < y < 5 m,
% 1,8 m < z < 5 m, D = 2(yz ax + xz ay ? 2xy az
% )/z2 C/m2
% . Avalie, numericamente, a
% carga dentro do volume usando dois métodos diferentes.

clc
clear all
close all

x1=2;
x2=3;
y1=2;
y2=3;
z1=2;
z2=3;
%k=1/(4*pi*8.854*10^(-12));  %Constante


passe = 1/500;
%Gerador do campo:
xl=[x1:passe:x2]; % variação da coordenada x onde está a carga 
yl=[y1:passe:y2];
zl=[z1:passe:z2];

V = ((passe)^3); %volume de cada segmento
Q = 0;

syms x y z 
D =[2*y/z 2*x/z -4*x*y/(z^2)]; %input('Insira o vetor na forma simbolica [fx(x,y,z) fy(x,y,z) fz(x,y,z)]...\n')
divD = diff(D(1),x) + diff(D(2),y) + diff(D(3),z);
pv=divD;

div = symfun(pv,[x y z]); %Cria Função Simbólica e substitui no Divergente
div2 = matlabFunction(div); %Pega funçao simbolica e passa pra numerica de novo(Acelera o processo)

for m = 1:length(xl)  % varre a coordenada x da carga
    disp(m)
    for n = 1:length(yl)  % varre a coordenada y da carga
        for o = 1:length(zl)

            Q = Q + (div2(xl(m),yl(n),zl(o))*V);  
            
        end
    end
end
disp(double(Q))

