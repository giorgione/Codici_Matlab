%ProvaWhitening.m
%
% Data una Distribuzione Gaussiana Bidimensionale con variabili
% statisticamente indipendenti, gli applico la trasformazione de whitening.
%
% Procedimento:
%
% 1) Costruisco F, la Distribuzione Gaussiana 2D
%
% 2) Calcolo la Matrice di Trasformazione Whitening
%
% 3) Calcolo Fw la Distribzione Trasformata
%
% 4) Disegno F per curve di livello: ellissi
%
% 5) Disegno Fw per curve di livello: circonferenze
clc;clear;close all;
syms x  y ;
X=[x;y];
M=[3;2];
%Matrice di Varianza Covarianza
S=[5^2 0;
   0 3^2];

P=GaussianaMulti(S,M,X);



%Applico la trasformazione whitening alla Distribuzione Gaussiana F: le sue
%curve di livello saranno circonferenze
[Media,Sigma,Aw]=WhiteningTransform(M,S);

[x y]=meshgrid(-10:.25:10,-10:.25:10);
[m,n]=size(x);
X=[reshape(x,1,m*n);reshape(y,1,m*n)];
Fw=zeros(m,n);

for i=1:m*n
    %valuto la distribuzione trasformata
    Fw(i)=GaussianaMulti(Sigma,Media,X(:,i));
    
    %valuto la distribuzione originale
    F(i)=GaussianaMulti(S,M,X(:,i));
end
Fw=reshape(Fw,m,n);
F=reshape(F,m,n);

xlabel('x');
ylabel('y');
contour(x,y,F)
hold on;
plot(M(1),M(2),'or','MarkerSize',3);

pause;
xlabel('x');
ylabel('y');
contour(x,y,Fw)

plot(Media(1),Media(2),'or','MarkerSize',3);


figure;
surf(x,y,F);title('Distribuzione Iniziale')

figure;
surf(x,y,Fw); title('Distribuzione trasformata mediante whitening')

