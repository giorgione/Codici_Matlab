% function F=GaussianaMulti_Fun(Sigma,M,X)
% Costruisce ( funzione Anonima da poter passare ad eventuali operazioni di
% Analisi) la Distribuzione Gaussiana Multivariata specificata dai
% parametri:
%
% Sigma : matrice di Varianza Covarianza
%
% M : Vettore delle Medie
%

function [fun]=GaussianaMulti_Fun(Sigma,M)
SigmaInv=inv(Sigma);
d=length(M);
%definisco la funzione anonima
fun= @(x) exp( (-.5)*( (x-M).'*SigmaInv*(x-M)) ) / (((2*pi).^(d/2))*(det(Sigma)^.5));