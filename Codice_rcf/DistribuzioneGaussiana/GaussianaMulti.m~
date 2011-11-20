% function F=GaussianaMulti(Sigma,M,X)
% Valuta la Distribuzione Gaussiana Multivariata specificata dai
% parametri:
%
% Sigma : matrice di Varianza Covarianza (le varianze sono elevate al quadrato)
%
% M : Vettore delle Medie
%
% X : Vettore delle coordinate del punto in cui si desidera valutare la Probabilit� di
% estrarre il campione X
function [F]=GaussianaMulti(Sigma,M,X)
SigmaInv=inv(Sigma);
%F=exp( (-.5)*( (X-M).'*SigmaInv*(X-M)) ) / (((2*pi).^(d/2))*(det(Sigma)^.5));
[d,N]=size(X);

Denum= ( (2*pi)^(d/2) )*(sqrt(det(Sigma)) );

for i=1:N
    F(i)=exp( -(( ( X(:,i)-M).'*SigmaInv*(X(:,i)-M)) )/2 );
end
F=F/Denum;