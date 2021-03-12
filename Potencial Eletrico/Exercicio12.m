%% Descobrir o E1 - Campo eletrico dentro do cilindro
%tem que ser por gauss pq é um fio infinito
 

%% Variaveis Dadas
clc
clear all
close all

%% variaveis do problema
e0 = 8.854*10^-12;
ps = 7.9*10^-12;

const=1/(4*pi*e0);  %Constante

%% Variaveis Criadas

passe = 0.0001;
limites =10;

%Gerador do campo:
xl= 0:passe:3.7; % variação da coordenada x onde está a carga 
yl= -7.1:passe:0; % variação da coordenada y onde está a carga 

dS = passe^2; %tamanho de cada segmento

%% Desenvolvimento
syms x z
V = (100*z);

V1 = symfun(V,z); %Cria Função Simbólica e substitui no Divergente
V2 = matlabFunction(V1);
E=0;
for i=0:1:10
    E = (V2(z+passe)-V2(z))./(passe);
end
E = -E*x/((x^2)+4);
D2 = E*e0;
D1 = symfun(D2,x); %Cria Função Simbólica e substitui no Divergente
D = matlabFunction(D1);
disp(D);

%Calcular carga interna
Q = 0;

for m = 1:length(xl)  % varre a coordenada x da carga
    for n = 1:length(yl)  % varre a coordenada y da carga
        Q = Q + double(D(xl(m))*dS);
    end       
end

disp(Q)