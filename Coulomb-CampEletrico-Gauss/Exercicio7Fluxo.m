%% Quest�o 7
% Na regi�o do espa�o livre que inclui o volume 1,8 m < x < 5 m, 1,8 m < y < 5 m,
% 1,8 m < z < 5 m, D = 2(yz ax + xz ay ? 2xy az)/z2 C/m2
% . Avalie, numericamente, a
% carga dentro do volume usando dois m�todos diferentes.

clc
clear all
close all

x1=1.8;
x2=5;


passo = 1/500;

%Gerador do campo:
xl=[x1:passo:x2]; % varia��o da coordenada onde est� a carga 
yl=[x1:passo:x2];
zl=[x1:passo:x2];

A = passo^2;
ss = 0; si = 0;

%Como � um paralelep�pedo, temos 6 Normais. 3 Eixos com coordenadas
%positivas e negativas.

syms x y z 
D =[2*y/z 2*x/z -4*x*y/(z^2)]; 

dens = symfun(D,[x y z]); %Cria Fun��o Simb�lica e substitui na Densidade
densN = matlabFunction(dens);   %Passa a fun��o de Simb�lica pra Num�rica


%Superf�cie Superior
dS = [0 0 1];
for i = 1:length(xl)  % varre a coordenada x da carga
    disp(i)
    for j = 1:length(yl)  % varre a coordenada y da carga
            ss = ss + A*double(dot((densN(xl(i),yl(j),x2)),dS)); 
    end
end
%Superf�cie Inferior
dS = [0 0 -1];
for i = 1:length(xl)  % varre a coordenada x da carga
    disp(i)
    for j = 1:length(yl)  % varre a coordenada y da carga
            si = si + A*double(dot((densN(xl(i),yl(j),x1)),dS)); 
    end
end


Qz = ss + si;
disp(double(Qz))
%RESPOSTA: 127.9341
