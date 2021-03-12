
% Na regi�o do espa�o livre que inclui o volume 1,8 m < x < 5 m, 1,8 m < y < 5 m,
% 1,8 m < z < 5 m, D = 2(yz ax + xz ay ? 2xy az
% )/z2 C/m2
% . Avalie, numericamente, a
% carga dentro do volume usando dois m�todos diferentes.

clc
clear all
close all

x1=1;
x2=2;
y1=1;
y2=2;
z1=1;
z2=2;
k=1/(4*pi*8.854*10^(-12));  %Constante


passe = 1/500;
%Gerador do campo:
xl=[x1:passe:x2]; % varia��o da coordenada x onde est� a carga 
yl=[y1:passe:y2];
zl=[z1:passe:z2];

V = ((passe)^3); %volume de cada segmento
Q = 0;

syms x y r 
D = [16*r*cos(2*x) 0 0];
divD = diff(D(1),x) + diff(D(2),y) + diff(D(3),z);
pv=divD;

div = symfun(pv,[x y r]); %Cria Fun��o Simb�lica e substitui no Divergente
div2 = matlabFunction(div); %Pega fun�ao simbolica e passa pra numerica de novo(Acelera o processo)

for m = 1:length(xl)  % varre a coordenada x da carga
    disp(m)
    for n = 1:length(yl)  % varre a coordenada y da carga
        for o = 1:length(zl)

            Q = Q + (div2(xl(m),yl(n),zl(o))*V);  
            
        end
    end
end

disp(double(Q))

