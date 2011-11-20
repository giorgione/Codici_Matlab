% function F=GaussianaMulti(Sigma,M,X)
% Valuta la Distribuzione Gaussiana Multivariata specificata dai
% parametri:
%
% Sigma : matrice di Varianza Covarianza
%
% M : Vettore colonna delle Medie
%
% X : Vettore delle coordinate dei punti in cui si desidera valutare la Probabilità di
% estrarre il campione X
function [F]=GaussianaMulti_Punti(Sigma,M,X)
SigmaInv=inv(Sigma);
[d,n]=size(X);
F=0;
for i=1:n
    F=F+log(exp( (-.5)*( (X(:,i)-M).'*SigmaInv*(X(:,i)-M)) ) / (((2*pi).^(d/2))*(det(Sigma)^.5)));
end