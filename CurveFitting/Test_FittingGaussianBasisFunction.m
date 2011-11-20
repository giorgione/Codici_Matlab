%Esempio di Curve Fitting con Gaussian Basis function
clc;clear;close all

x=linspace(0,1,50);
%curva da stimare
y=sin(2*pi*x);

%Valori osservati con Rumore additivo
y1=randperm(50);

%Numero punti osservati
n=10;
rumore=randn(1,n);
t=y(y1(1:n))+rumore;
xt=x(y1(1:n));

plot(x,y,'r');hold on
plot(xt,t,'ob');


%Approssimazione attraverso Polinomio di grado 3
basi=@(x,u,s)([ (x-u(1))/(2*s^2) , (x-u(2))/(2*s^2) ,(x-u(3))/(2*s^2)]);
%calcolo i coefficienti del polinomio
[w]=LinearBasis(t',xt,basi)


% Valuto il polinomio su tutto il dominio di rappresentazione
for i=1:50
    y(i)=sum(basi(x(i))*w);
    
end
plot(x,y,'g')

% La distribuzione predittiva sarà
% N(t| y(x,w), b)
%
w_ml=w;

%Calcolo la varianza della Gaussiana
s=0;
for i=1:n
    s=s+( t(i)-basi(xt(i))*w_ml )^2;
end
s=s/n;

s=sqrt(1/s);
%Valuto la distribuzione predittiva su t che mi restituisce la probabilità che t
% assuma un determinato valore in corrispondenza di xt(i)
for i=1:10
    Media=basi(xt(i))*w_ml;
    
    % Calcolo la distribuzibe su
    T=linspace(Media-4*s,Media+4*s,50);
    
    N= exp(-( (Media-T).^2 )./(2*s^2) ) ./sqrt(2*pi*s^2);
    
       
    plot(xt(i)+N,T,'m');
end