% function DisegnaEllisse(Xo,rx)
% Disegna l' Ellisse di centro Xo=[x y] e assi rx e ry
% Parametri in:
% 
% - Xo: centro dell Ellisse
%
% - rx: Asse x dell Ellisse
%
% - ry: Raggio dell Ellisse
%
% - n: Numero di Punti in cui si desidera valutare l Ellisse
% 
% Parametri out:
%
% - x: n ascisse sull Ellisse
%
% - y: n ordinate sull Ellisse
function [x,y]=DisegnaEllisse(Xo,rx,ry,n,disegna)
t=linspace(0,2*pi,n);
x=Xo(1)+rx*cos(t);
y=Xo(2)+ry*sin(t);
if disegna==1 
    plot(x,y,'r')
end