%% Questão 7
% Na região do espaço livre que inclui o volume 1,8 m < x < 5 m, 1,8 m < y < 5 m,
% 1,8 m < z < 5 m, D = 2(yz ax + xz ay ? 2xy az)/z2 C/m2
% . Avalie, numericamente, a
% carga dentro do volume usando dois métodos diferentes.

clc
clear all
close all

x1=1.8;
x2=5;


passo = 1/500;

%Gerador do campo:
xl=[x1:passo:x2]; % variação da coordenada onde está a carga 
yl=[x1:passo:x2];
zl=[x1:passo:x2];

A = passo^2;
ss = 0; si = 0;

%Como é um paralelepípedo, temos 6 Normais. 3 Eixos com coordenadas
%positivas e negativas.

syms x y z 
D =[2*y/z 2*x/z -4*x*y/(z^2)]; 

dens = symfun(D,[x y z]); %Cria Função Simbólica e substitui na Densidade
densN = matlabFunction(dens);   %Passa a função de Simbólica pra Numérica


%Superfície Superior
dS = [0 0 1];
for i = 1:length(xl)  % varre a coordenada x da carga
    disp(i)
    for j = 1:length(yl)  % varre a coordenada y da carga
            ss = ss + A*double(dot((densN(xl(i),yl(j),x2)),dS)); 
    end
end
%Superfície Inferior
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
