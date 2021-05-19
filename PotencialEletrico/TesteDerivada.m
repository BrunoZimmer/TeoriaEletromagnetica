clc
clear all
close all

y = @(x) x.^4;
value = 2;
dx = 0.01;

for i=0:dx:2
    yl = (y(i+dx)-y(i))./(dx);
end
disp(yl);

syms x
%y4 = x.^4;

y4= (100*x)/((x^2)+4);

y2 = symfun(y4,x); %Cria Função Simbólica e substitui no Divergente
y2 = matlabFunction(y2);
y3=0;
for i=0:dx:2
    y3 = (y2(x+dx)-y2(x))./(dx);
end
disp(y3);
disp(double(subs(y3,x,2)));