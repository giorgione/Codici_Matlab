% function F=GaussianaMulti_Fun_t(t)
% Costruisce ( funzione Anonima da poter passare ad eventuali operazioni di
% Analisi) la Distribuzione Gaussiana Multivariata specificata dai
% parametri:
%
% Sigma : matrice di Varianza Covarianza
%
% M : Vettore delle Medie
%
% X : Vettore dei punti in cui si desidera valutare la Probabilità di
% estrarre il campione costituito dai punti in X
function [fun]=GaussianaMulti_Fun_t(t)
Sigma=diag(t(3:4));
SigmaInv=inv(Sigma);
d=2;
%definisco la funzione anonima
fun= @(x) exp( (-.5)*( (x-t(1:2)).'*SigmaInv*(x-t(1:2))) ) / (((2*pi).^(d/2))*(det(Sigma)^.5));