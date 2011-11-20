function DisegnaPiano(a,b)
%Disegna lo span(a,b): spazio generato dai Vettori a e b
syms x y 
P=a*x+b*y;
ezsurf(P)
hold on;
compas3d(a(1),a(2),a(3),'r')
compas3d(b(1),b(2),b(3),'r')
 
 
 