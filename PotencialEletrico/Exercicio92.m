%% Descobrir o E1 - Campo eletrico dentro do cilindro
%tem que ser por gauss pq é um fio infinito
 

%% Variaveis Dadas
clc
clear all
close all

%variaveis do problema
e0 = 8.854*10^-12;
a = 2.7*10^-3;
b = 7.3*10^-3;

pv = 8.1*10^-12;
ps = -pv*((a^2)/(2*b));

k=1/(4*pi*e0);  %Constante

%% Variaveis Criadas

passe = a/5;
limites = b*2;
%Onde o campo sera medido:
x=[-limites:passe:limites]; %vetor na coordenada x onde será calculado E
y=[-limites:passe:limites]; %vetor na coordenada z onde será calculado E
z=[-limites:passe:limites]; %vetor na coordenada z onde será calculado E

%Gerador do campo interno (a):
xl=[-a:passe:a]; % variação da coordenada x onde está a carga 
yl=[-a:passe:a];
zl=[-limites:passe:limites];

%Gerador do campo casca externa:
xll=[-b:passe:b]; % variação da coordenada x onde está a carga 
yll=[-b:passe:b];
zll=[-limites:passe:limites];

dV = passe^3; %volume de cada segmento
dS = passe^2; %área de cada segmento

%inicializa o campo elétrico:
Va = zeros (length(x),length(y)); 
Vb = zeros (length(x),length(y)); 

%% Desenvolvimento

%Campo dentro de A
for i = 1:length(x)% varre a coordenada x onde E será calculado
    disp(i)
    for j = 1: length(y)  % varre a coordenada y onde E será calculado

        for m = 1:length(xl)  % varre a coordenada x da carga
            for n = 1:length(yl)  % varre a coordenada y da carga

                r = [x(i),y(j),0]; %vetor posição apontando para onde estamos calculando E
                rl= [xl(m),yl(n),0];% vetor posição apontando para um segmento do disco
                %disp((r-rl)*(r-rl)')
                if ((r-rl)*(r-rl)'>passe/1000) && (sqrt(xl(m)^2+yl(n)^2) <= a)
                    Va(i,j) = Va(i,j)  + k*((pv*dV)/sqrt((r-rl)*(r-rl)'))';
                end  

%                 if (sqrt(x(i)^2+y(j)^2) > b)
%                     Va(i,j) = 0;
%                 end
            end
        end
            
    end
end

%Campo entre A e B
for i = 1:length(x)% varre a coordenada x onde E será calculado
    disp(i)
    for j = 1: length(y)  % varre a coordenada y onde E será calculado

        for m = 1:length(xll)  % varre a coordenada x da carga
            for n = 1:length(yll)  % varre a coordenada y da carga

                r = [x(i),y(j),0]; %vetor posição apontando para onde estamos calculando E
                rll= [xll(m),yll(n),0];% vetor posição apontando para um segmento do disco

                if ((r-rl)*(r-rl)'>passe/1000)
                    if (sqrt(xll(m)^2+yll(n)^2) >= 0.96*b) && (sqrt(xll(m)^2+yll(n)^2) <= 1.04*b)
                        Vb(i,j) = Vb(i,j)  + k*((ps*dS)/sqrt((r-rll)*(r-rll)'))';
                    end
                end  

%                 if (sqrt(x(i)^2+y(j)^2) > b)
%                     Vb(i,j) = 0;
%                 end

            end
        end
    end
end
%% Grafico
% xd = linspace(-0.1,0.1);
% yd = linspace(-0.1,0.1);
% [X,Y] = meshgrid(xd,yd);

V = Va + Vb;

figure(1)
title('Va', 'Color', 'b')
surf(x,y,Va);
xlabel('x')
ylabel('y')
zlabel('z')
%axis([-0.1 0.1 -0.1 0.1 0.1 0.000004])
grid on
colormap(jet(20))
colorbar
hold on

figure(2)
title('Vb', 'Color', 'b')
surf(x,y,Vb);
xlabel('x')
ylabel('y')
zlabel('z')
%axis([-0.1 0.1 -0.1 0.1 0.1 0.000004])
grid on
colormap(jet(20))
colorbar

figure(3)
title('Va * Vb', 'Color', 'b')
surf(x,y,V);
xlabel('x')
ylabel('y')
zlabel('z')
%axis([-0.1 0.1 -0.1 0.1 0.1 0.000004])
grid on
colormap(jet(20))
colorbar


max = int64(length(xl));
resp = 0;

V0 =  V(max/2,max/2);
disp(0-V0)
disp(resp)
%0,62