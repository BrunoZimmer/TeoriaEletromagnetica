%% Questao 2
% Uma densidade volum�trica de carga uniforme de 0,2 ?C/m3 est� presente em
% uma casca esf�rica que se estende de r = 3 cm a r = 7,4 cm. Se pv = 0 em qualquer
% outra regi�o, calcule numericamente a carga total presente na casca.

clc
clear all
close all

fator = 10000;
rmax = 7.4/100; %m
disp(rmax);
rmin = 3/100; %m
disp(rmin);
passo = 1/fator;
densidade = 0.2*(10^-6); %uC/m3
dq = (passo^3)*densidade;
tamanho = int16(rmax*fator + 1);
carga = 0;

xl=[0:1:tamanho];
yl=[0:1:tamanho];
zl=[0:1:tamanho];

for x = 1:length(xl)  % varre o �ngulo da casca
    x
    for y = 1: length(yl)  % varre o �ngulo da casca
        for z = 1:length(zl)% varre o raio onde se encontra a carga
            raio = sqrt(x^2+y^2+z^2)*passo;
            if (raio <= rmax) && (raio >= rmin)
                carga = carga + dq;
            end
        end
    end
end

disp(double(carga)*8);


