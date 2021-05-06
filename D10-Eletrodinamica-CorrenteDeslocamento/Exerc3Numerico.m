% Uma espira filamentar quadrada de lado a está inicialmente posicionada
% no plano xy, com seu centro alinhado com o eixo z, numa região de campo
% uniforme B(y,t) = B0 az. A espira começa a girar com velocidade angular
% ? sobre o eixo x. Determine numericamente o fluxo magnético através da
% espira e a FEM induzida nessa espira ao longo do tempo. Trace os gráficos
% das duas quantidades em função do tempo.

clear all
clc
close all

%% Variaveis dadas

u0=4*pi*10^(-7); % Permeabilidade magnÃ©tica do espaÃ§o livre (em H/m)
B0 = 5.3;
w= 5;
a= 2;

%% Espaços de calculos

passo=a/5; % passo

dx=passo;
dy=passo;
dz=passo;

x= -2*a:dx:2*a; %variacao da coordenada x
y= -2*a:dy:2*a; %variacao da coordenada y
z= -2*a:dz:2*a; %variacao da coordenada z

xmedio = ceil(length(x)/2);
ymedio = ceil(length(y)/2);
zmedio = ceil(length(z)/2);

xl= -a/2:dx:a/2; %variacao da coordenada y onde está a placa
yl= -a/2:dy:a/2; %variacao da coordenada z onde está a placa

pl = 0:dx/2:a/2;

t = 0:pi/10:4*pi;

%% Descobrindo valor do fluxo magnético


% Determine a densidade de fluxo magnético

%zerando a densidade de fluxo magnético
% Flx(:,:,:,:,:) = zeros(3,length(x),length(y),length(z), length(t));
Flx(:,:) = zeros(3, length(t));

for p = 1:length(t)% varre a coordenada x onde H será¡ calculado
%     for i = 1:length(x)% varre a coordenada x onde H será¡ calculado
%         disp(i)
%         for j = 1: length(y)  % varre a coordenada y onde H será calculado
%             for k = 1:length(z) % varre a coordenada z onde H será calculado
%                 t = 0; %definir nesse tempo o calculo pra depois tentar fazer girar
                for m = 1:length(xl) 
                    for n = 1:length(yl)
%                 for o = 1:length(pl)

%                         r = [x(i),y(j),z(k)];  %vetor posiÃ§Ã£o do ponto de campo (onde calculamos B)
                        B = [0,0,B0];
%                         rhol1 = sqrt(xl(m)^2+yl(n)^2);
%                         [rl1(1), rl1(2), rl1(3)] = pol2cart(t(p), pl(o), 0); % convertemos rl para coord retangulares
                        an = [0, -1*sin(w*t(p)), cos(w*t(p))];

                        dS = dx*dy; %componente de espaço bidimensional cart.
                        % Utilizamos a condiÃ§Ã£o abaixo para evitar a divisÃ£o por zero, no caso em que r e rl sÃ£o iguais.
%                         if ((r-rl1)*(r-rl1)' > (dz/4)^2)
                            dSn = dot(an,dS);
                            dFlx_ijk1 = dot(B,dSn); % essa  a contribuicao dB devida a esse dL

                            Flx(1,p) = Flx(1,p) + dFlx_ijk1(1);
                            Flx(2,p) = Flx(2,p) + dFlx_ijk1(2);
                            Flx(3,p) = Flx(3,p) + dFlx_ijk1(3);
%                         end

                    end
                end
%             end
%         end
%     end
    
end

%iniciando o frame
% init_getframe = struct('cdata',[],'colormap',[]);
% frames = repmat(init_getframe, length(t), 1 );
% frames(1) = getframe;
% h = quiver(Y, Z, squeeze(B(2,xmedio,:,:,:))', squeeze(B(3,xmedio,:,:,:))');

% for i = 2:length(t)
%     set(h,'XData',t(1:i));
%     set(h,'YData',B(1:i));
%     drawnow
%     frames(i) = getframe;


figure(1)
[Y,Z] = meshgrid(y,z);

% end
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';

frames = zeros(length(t)-1);

for i = 2:length(t)
    h = quiver(Y, Z, squeeze(B(2,xmedio,:,:,i))', squeeze(B(3,xmedio,:,:,i))');
    drawnow 
    % Capture the plot as an image 
    frame = getframe;
%     frames(i) = frame;
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    
    if i == 2 
        imwrite(imind,cm,filename,'gif', 'Loopcount',0); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end 
end



