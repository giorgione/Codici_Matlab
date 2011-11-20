% DistribuzioneGaussianaMulti_sc.m
% 1) Costruisco una Distribuzione Gaussiana Bivariata a variabili statisticamente
% Indipendenti
%
% 2) Disegno la Distribuzione Gaussiana
clc;clear;
syms x y Ux Uy Ox Oy Oxy
X=[x;y];
M=[Ux;Uy];

%Matrice di Varianza Covarianza
S=[Ox^2 Oxy;
   Oxy Oy^2];

Sinv=inv(S);
d=2;

Px=1/((2*pi)^(d/2)*det(S)^.5)*exp(-1/2*(X-M).'*Sinv*(X-M));

%Considero x1 e x2 statisticamente indipendenti: Oxy=0
Px1=subs(Px,Oxy,0);
Px1=simplify(Px1);
pretty(Px1);

% Distribuzione Bivariata con:
%
% medie :Ux=2  Uy=1
%
% Varianza-covarianza  Oy=.5 Ox=1 Oxy=Oyx=1
%
% x e y statisticamente indipendenti --> Oxy=0

Px1=subs(Px1,{Ux,Ox,Uy,Oy},{2,1,1,.5});
ezsurf(Px1)