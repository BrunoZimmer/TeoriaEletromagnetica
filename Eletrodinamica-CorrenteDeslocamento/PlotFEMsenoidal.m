clear all
clc
close all

%Variaveis dadads
B0 = 5.3;
w= 5;
a= 0.02;
t = 0:pi/100:2*pi;
y = (a^2)*B0*w*sin(w*t); %y = FEM
Sizet = length(t);

% Plot
figure(1)
h = plot(t(1),y(1));
xlim([t(1) t(end)])
ylim([min(y) max(y)])
xlabel('tempo')
ylabel('FEM')
legend('Magnitude da FEM induzida (âz)')
init_getframe = struct('cdata',[],'colormap',[]);
frames = repmat(init_getframe, Sizet, 1 );
frames(1) = getframe;

% Projetando frames
for i = 2:Sizet
    set(h,'XData',t(1:i));
    set(h,'YData',y(1:i));
    drawnow
    frames(i) = getframe;
end
% Play movie
movie(frames)