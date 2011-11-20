%Fitting della curva t=f(x)+e
%
% t: target values osservati
%
% x: input values
%
% basi: linear basis function ( basi=@(x)([1 x x^2 x^3]))
%
% w: coefficienti del Modello Lineare caratterizzato dalle basi

function [w,DesignMatrix]=LinearBasis(t,x,basi)
n=length(x);
To=basi(x(1));
%Determino il numero delle Basi
[N M]=size(To);

%Calcolo la Matrice PseudoInversa (n x M)
DesignMatrix=zeros(n,M);
for i=1:n
    DesignMatrix(i,:)=basi(x(i));
end
It=DesignMatrix';

%Calcolo la stima ML dei coefficienti
P=inv(It*DesignMatrix)*It;
w=P*t;
%w=w_ml;

%b_ml=sum((I*w-t).^2)/n;
%b=b_ml;


