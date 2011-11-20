% DisegnaGaussiana1_sc.m
% 1) Costruisco una Distribuzione Gaussiana Monovariata Normale
%
% 2) Disegno la Distribuzione Gaussiana Normale

clc;clear;close all

Ux=1;
Ox=2;

M=[Ux];

x=linspace(-8,8,150);

%Matrice di Varianza Covarianza
S=[Ox^2];

[m,n]=size(x);

F=zeros(1,m*n);


F=GaussianaMulti(S,M,x);

plot(x,F);
hold on
plot(Ux,0,'bo','MarkerFaceColor','r');

Area_Tot=sum(F)