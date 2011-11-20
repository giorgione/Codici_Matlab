%ProvaWhitening_1_3dim.m
%
% Data una Distribuzione Gaussiana Tredimensionale con variabili
% statisticamente indipendenti, gli applico la trasformazione de whitening.
%
% Procedimento:
%
% 1) Costruisco F, la Distribuzione Gaussiana 3D
%
% 2) Calcolo la Matrice di Trasformazione Whitening
%
% 3) Disegno Fw per curva di livello caratterizzata dalla distanza di Malanobhis
%    costante: sfera
%
clc;clear;close all;
syms x  y z;
X=[x;y;z];
M=[1;4;8];
%Matrice di Varianza Covarianza Sigma=I
S=[1^2 0  0;
   0 1^2  0;
   0  0  1^2];

P=GaussianaMulti(S,M,X);

%Applico la trasformazione whitening alla Distribuzione Gaussiana F
[Media,Sigma,Aw]=WhiteningTransform(M,S);


hold on;
plot3(M(1),M(2),M(3),'or','MarkerSize',3);
plot3(Media(1),Media(2),Media(3),'or','MarkerSize',3);
syms x y;
Pw=GaussianaMulti(Sigma,Media,X);

%Calcolo la Sfera caratterizzata da tutti i feature-Vector che hanno
%distanza r=3 dal Vettore delle Medie
% Poich� la varianza di x,y,z � ro=1 ottengo che tutti il 99% dei feature
% vector cadono all' interno della sfera di raggio 3*ro calcolata tramite
% l' equazione della distanza di Mahalanobis
r=3;
Mal=r^2-(X-Media).'*inv(Sigma)*(X-Media);
Sol=solve(Mal,'z');

z=3*cplxgrid(20);
X=M(1)+real(z);
Y=M(2)+imag(z);

F1=8+(-8-X.^2+2*X-Y.^2+8*Y).^(1/2)
F2= 8-(-8-X.^2+2*X-Y.^2+8*Y).^(1/2)
figure;
surf(X,Y,F1);hold on;
surf(X,Y,F2)
axis tight

