% EsempioMinimiQuadrati.m
% Problema ai minimi quadrati su 4 Punti in R2.
% Procedura  one-shot applicata sull ' intero Set di Dati

clc;clear;close all
%Considero i 2 punti di R1
P1=[1 2];
P2=[2 0];

%Considero i 2 punti di R2
P3=[3 1];
P4=[2 3];
Punti=[P1; P2;P3;P4];
plot(Punti(1:2,1), Punti(1:2,2),'or','MarkerSize',5,'MarkerFaceColor','r')
hold on;
plot(Punti(3:4,1), Punti(3:4,2),'ob','MarkerSize',5,'MarkerFaceColor','b')
%Creo la Matrice Y per l' applicazione dei minimi quadrati

Y=[1   P1;
   1   P2;
   -1 -P3;
   -1 -P4];

%Considero il vettore  b a coefficienti positivi
b=[1;1;1;1];

%Calcolo la pseudo inversa di Y
Ypsi=inv(Y.'*Y)*Y.';

a=Ypsi*b;

syms x1 x2
F=a.'*[1;x1;x2];
S=solve(F,'x2');
ezplot(S)
axis equal
