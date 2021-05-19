% Uma espira filamentar quadrada de lado a está posicionada no plano xy, 
% com seu centro alinhado com o eixo z. A espira está parada nesta posição.
% Determine a magnitude da FEM induzida que resulta de uma densidade de 
% fluxo magnético dado por B(y,t) = B0 cos( ?t - ?y) az. Afim de avaliar 
% a sua resposta considere B0 = 8,1 Wb/m2,? = 54,9 rad/s, ? = 2,4 rad/m, 
% a = 0,01 m, t = 2,8 s.


clear all
clc
close all

%% Variaveis dadas

u0=4*pi*10^(-7); % Permeabilidade magnÃ©tica do espaÃ§o livre (em H/m)
a = 1;
B0 = 8.1;
w = 54.9;
beta = 2.4;
tmp = 2.8;


%% Espaços de calculos

passo=a/5; % passo

dx=passo;
dy=passo;
dz=passo;

x= -2*a:dx:2*a; %variacao da coordenada x
y= -2*a:dy:2*a; %variacao da coordenada y
z= -2*a:dz:2*a; %variacao da coordenada z

xl= -a/2:dx:a/2;
yl= -a/2:dz:a/2;

xmedio = ceil(length(x)/2);
ymedio = ceil(length(y)/2);
zmedio = ceil(length(z)/2);


f=w/(2*pi);
T=1/f;

dt = T/100;
t = 0:dt:2*tmp;

%% Descobrindo valor do campo B

%zerando a densidade de fluxo magnético
B(:,:,:,:,:) = zeros(3, length(x), length(y), length(z), length(t));

for p = 1:length(t)
    for i = 1:length(x)
        disp(i)
        for j = 1:length(y)  
            for k = 1:length(z) 
                
                B(3,i,j,k,p) = B0*cos((w*t(p))-(beta*y(j)));
                
            end
        end
    end
    
end

%% Plot animado do campo B

figure(1)
    % PARA GERAR IMAGEM GIF
% filename = 'testAnimated.gif';
% 
% frames = zeros(length(t)-1);

for i = 2:length(t)
    [Y,Z] = meshgrid(y,z);
    h = quiver(Y, Z, squeeze(B(2,xmedio,:,:,i))', squeeze(B(3,xmedio,:,:,i))');
    drawnow 
    pause(0.1) 
    % PARA GERAR IMAGEM GIF
%     % Capture the plot as an image 
%     frame = getframe;
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256); 
%     
%     % Write to the GIF File 
%     if i == 2 
%         imwrite(imind,cm,filename,'gif', 'Loopcount',0); 
%     else 
%         imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%     end 
end

%% Descobrindo o fluxo

fluxo = zeros(1, length(t));
for i = 1:length(xl)
    disp(i)
    for j = 1:length(yl)
        for p = 1:length(t)
            
            dS=dx*dy;
            n = [0,0,1];
            
            dSn=n*dS;
            
            fluxo(1,p) = fluxo(1,p) + dot(B(:,i,j,zmedio,p),dSn);
        end
    end
end


%% Plot Flx x tmp
figure(2);
plot(t,fluxo);
xlabel('Tempo')
ylabel('W')
title('Flx x tmp');

%% Descobrindo a FEM
FEM = zeros(1,length(t));

for i = 2:length(t)-1% varre a coordenada x onde E será calculado
    disp(i);
    FEM(1,i) = (-1 * (fluxo(1,i+1)-fluxo(1,i-1))) / (2*dt);
end

%% Plot FEM x tmp
figure(3);
plot(t,FEM);
xlabel('Tempo')
ylabel('V')
title('FEM x tmp');

disp(FEM(length(t)));
disp(FEM(length(t)));

