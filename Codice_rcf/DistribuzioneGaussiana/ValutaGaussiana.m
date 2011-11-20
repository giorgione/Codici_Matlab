%function F=ValutaGaussiana(S,M,a,b,n)
% Valuta la Gaussiana specificata da S ed M nell' intervallo [a,b]
%
% parametri in:
%
%
function [F,x,y]=ValutaGaussiana(S,M,a,b,n)

t=linspace(a,b,n);
[x y]=meshgrid(t);

[m,n]=size(x);
X=[reshape(x,1,m*n);reshape(y,1,m*n)];
F=zeros(m,n);

for i=1:m*n
    F(i)=GaussianaMulti(S,M,X(:,i));
end
F=reshape(F,m,n);