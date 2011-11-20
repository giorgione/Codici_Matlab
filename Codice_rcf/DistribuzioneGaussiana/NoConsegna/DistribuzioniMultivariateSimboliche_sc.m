% DistribuzioneMultivariateSimboliche_sc.m
% Calcolo Simbolico Distribuzione Gaussiana Multivariata per variabili 
% statisticamente Indipendenti

clc;clear;
syms x y z Ux Uy Uz Ox Oy Oz Oxy Oxz Oyz
X=[x;y;z];
M=[Ux;Uy;Uz];
%Matrice di Varianza Covarianza
S=[Ox^2 Oxy Oxz;
   Oxy Oy^2 Oyz;
   Oxz Oyz Oz^2];

Sinv=inv(S);
d=length(X);

P=GaussianaMulti(S,M,X);

%Considero x,y,z statisticamente indipendenti: Oxy=0
P=subs(P,{Ux,Uy,Uz,Ox,Oy,Oz,Oxy,Oxz,Oyz},{0,0,0,1,1,1,0,0,0});
P=simplify(P);




