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
clc;clear;close all

Ux=2;
Ox=1;

Uy=1;
Oy=.5;

Oxy=0.2;

M=[Ux;Uy];

[x y]=meshgrid(-2:.125:8,-2:.125:8);

%Matrice di Varianza Covarianza
S=[Ox^2 Oxy;Oxy Oy^2];
d=2;
[m,n]=size(x);
X=[reshape(x,1,m*n);reshape(y,1,m*n)];
F=zeros(m,n);

for i=1:m*n
    F(i)=GaussianaMulti(S,M,X(:,i));
end
F=reshape(F,m,n);

%Calcola la probabilità condizionata P(x|y=0)
hold on;

mesh(x,y,F);
xlabel('x');
ylabel('y');

figure;
contour(x,y,F);
xlabel('x');
ylabel('y');

%Calcolo gli Autovalori ed Autovalori della matrice Sigma
[X,lamda]=eig(S);
hold on;
a1=sqrt(lamda(1,1));
a2=sqrt(lamda(2,2));

line([X(1,1)+Ux; -X(1,1)+Ux],[X(2,1)+Uy;-X(2,1)+Uy])
line([X(1,2)+Ux; -X(1,2)+Ux],[X(2,2)+Uy;-X(2,2)+Uy])
