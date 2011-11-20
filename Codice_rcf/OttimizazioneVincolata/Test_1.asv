%Test sul Problema do Minimizazione Vicnocolata
clear;close all;clc
syms x y;
f1=figure;
F=rosenbrock([x;y]);
ezsurf(F,[-50 50 -50 50]);
hold on;
ezcontour(F,[-50 50 -50 50]);
l=ezplot('x^2+y^2-1')
set(l,'Color','g','LineWidth',2);

options = optimset('Display','iter');
[x,fval] = fmincon(@rosenbrock,[0 0],[],[],[],[],[],[],@unitdisk,options);

plot3(x(1),x(2),fval,'or','MarkerSize',5,'MarkerFaceColor','r')

