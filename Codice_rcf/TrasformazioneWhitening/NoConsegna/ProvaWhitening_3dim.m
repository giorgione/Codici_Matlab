% Trasformo attraverso la trasformata Whitening una Distribuzione Normale a
% 3-dimensioni ()
%
clc;clear;close all;
syms x  y z;
X=[x;y;z];
M=[0;0;0];
%Matrice di Varianza Covarianza
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


