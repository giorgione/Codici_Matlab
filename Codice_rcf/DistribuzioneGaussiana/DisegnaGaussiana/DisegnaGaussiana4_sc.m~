% DisegnaGaussiana3_sc.m
% 1) Costruisco una Distribuzione Gaussiana Bivariata a variabili
%    correlate.
%
% 2) Disegno la Distribuzione Gaussiana Bivariata
%
% 3) Calcolo gli Auotovalori ed Autovalori della Matrice di Varianza
%    Covarianza e verifico come essi rappresentino rispettivamente la
%    lunghezza e direzione degli assi delle curve di livello della Guassiana
%
%
clc;clear;close all

syms Ux Uy Ox Oy Oxy x y p real;

%p=Oxy/(Ox*Oy);
M=[Ux;Uy];

%Matrice di Varianza Covarianza
S=[Ox^2      p*(Ox*Oy); 
   p*(Ox*Oy)     Oy^2 ];

d=2;

F=GaussianaMulti(S,M,[x;y])

C=log(.4*(2*pi*Ox*Oy*(1-p^2)^.5));
%Equazione Quadrica
Cl=((x-Ux)/Ox)^2+((y-Uy)/Oy)^2-2*p*((x-Ux)/Ox)*((y-Uy)/Oy)-C;

k=menu('Seleziona il Coefficiente di Correlazione:','1) 0','2) 0.99','3) -0.99');
switch k
    case 1
        P=0;
    case 2
        P=.99;
    case 3
        P=-.99;
end
Cl=subs(Cl,{Ux,Uy,Ox,Oy,p},{3,2,2,2,P})
ezplot(Cl,[-20 20])

S=double(subs(S,{Ux,Uy,Ox,Oy,p},{3,2,2,2,P}));
%Calcolo gli Autovalori ed Autovalori della matrice Sigma
[X,lamda]=eig(S);
hold on;
a1=lamda(1,1).^.5;
a2=lamda(2,2).^.5;

Ux=double(subs(Ux,3));
Uy=double(subs(Uy,2));
h1=line([a1*X(1,1)+Ux; -a1*X(1,1)+Ux],[a1*X(2,1)+Uy;-a1*X(2,1)+Uy]);
h2=line([a2*X(1,2)+Ux; -a2*X(1,2)+Ux],[a2*X(2,2)+Uy;-a2*X(2,2)+Uy]);

set(h1,'Color','r');
set(h2,'Color','r');

% F=GaussianaMulti(S,[3;2],[x;y])
% figure;
% ezsurf(F)