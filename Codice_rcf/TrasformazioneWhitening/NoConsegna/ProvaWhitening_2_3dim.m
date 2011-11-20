%ProvaWhitening_2_3dim.m
%
% Trasformo attraverso la trasformata Whitening una Distribuzione Gaussiana a
% 3-dimensioni traslata nel punto M(1,4,8) -> vettore delle medie 
% con Matrice di varianza covarianza Arbitraria (Ellissoide)
%
% Procedimento:
%
% 1) Costruisco F, la Distribuzione Gaussiana 3D
%
% 2) Calcolo la Matrice di Trasformazione Whitening e calcolo la
%    distribuzione trasformata Fw
%
% 3) Disegno F per curva di livello caratterizzata dalla distanza di Malanobhis
%    costante: ellissoide
%
% 3) Disegno Fw per curva di livello caratterizzata dalla distanza di Malanobhis
%    costante: sfera
%

clc;clear;close all;
syms x  y z;
V=[x;y;z];
M=[1;4;8];
rox=1.5;
roy=1;
roz=2;
%Matrice di Varianza Covarianza Sigma=I
S=[rox^2 0  0;
   0 roy^2  0;
   0  0  roz^2];

P=GaussianaMulti(S,M,V);

%Applico la trasformazione whitening alla Distribuzione Gaussiana con:
% M -> vettore delle Medie
% S -> Matrice di Varianza Covarianza
[Media,Sigma,Aw]=WhiteningTransform(M,S);
%La trasformata altera il Vettore delle Medie e la Matrice Sigma

figure;
hold on; grid on;
plot3(M(1),M(2),M(3),'or','MarkerSize',5,'MarkerFaceColor','r'); 
plot3(Media(1),Media(2),Media(3),'ob','MarkerSize',5,'MarkerFaceColor','b');
xlabel('x');
ylabel('y');
zlabel('z');

%Calcola l' iperEllissoide dei punti equidistanti di 1 dal Vettore delle Medie
%dis=max(max(S(1:2,1:2).^.5));
r=3;
Mal1=r^2-(V-M).'*inv(S)*(V-M);
Sol1=solve(Mal1,'z');

%Disegna l' ellissoide
[X,Y]=GrigliaEllittica(3*rox,3*roy,19);
X=X+M(1);
Y=Y+M(2);


F1=8+2/3*(-67-4*X.^2+8*X-9*Y.^2+72*Y).^(1/2);
F2=8-2/3*(-67-4*X.^2+8*X-9*Y.^2+72*Y).^(1/2);


figure;
surf(X,Y,F1);hold on;
surf(X,Y,F2)

pause;
%Calcolo la Sfera caratterizzata da tutti i feature-Vector che hanno
%distanza r=3 dal Vettore delle Medie
% Poichè la varianza di x,y,z è ro=1 ottengo che tutti il 99% dei feature
% vector cadono all' interno della sfera di raggio 3*ro calcolata tramite
% l' equazione della distanza di Mahalanobis

%Calcola la Sfera dei punti equidistanti dal Vettore delle Medie
%Trasformate mediante Whitening
r=3;
Mal2=r^2-(V-Media).'*inv(Sigma)*(V-Media);
Sol2=solve(Mal2,'z');

z=3*cplxgrid(20);
X=Media(1)+real(z);
Y=Media(2)+imag(z);


F1= 4+1/3*(-67-9*X.^2-72*X-9*Y.^2-12*Y).^(1/2);
F2= 4-1/3*(-67-9*X.^2-72*X-9*Y.^2-12*Y).^(1/2);

surf(X,Y,F1);hold on;
surf(X,Y,F2)
