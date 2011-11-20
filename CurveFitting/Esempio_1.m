% Esempio di Fitting attraverso minimizzazione della funzione obbiettivo
% SUM OF SQUARED ERROR SQR
%         
%                                                              M
% E(w)=    sqrt( 2*E(w)/N)
%   rms    
%      
%
% UTilizzo polinomi di grado M=2 - 3
clc;clear;close all

%%% Genero la funzione da Predire
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

syms w1 w2 w3 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Polinomio di grado M=2
W=[w1;w2; w3];

%Calcolo la Funzione obbiettivo: E(w)=(w1+w2*x(i)+w3*x(i)^2 -t(i))^2
Ew=0;
for i=1:n
    X=[1;xt(i);xt(i)^2];
    Ew=Ew+(W.'*X -t(i))^2;
end
Ew=Ew/2;

Erms=sqrt(2*Ew/n);

%Calcolo la derivata dEw/dw di Ew rispetto ai parametri w
DErms=jacobian(Erms,[w1 w2 w3])

%Calcolo la soluzione analitica risolvendo dEw/dw=0
W_analitica=solve(DErms);
W1=[double(W_analitica.w1);double(W_analitica.w2);double(W_analitica.w3);]

%Valuto la soluzione su tutti i punti x
for i=1:50
    X=[1;x(i);x(i)^2];
    y_stimata(i)=W1.'*X;
end

plot(x,y_stimata,'m');hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Polinomio di grado M=3
syms w4 

W=[w1;w2; w3;w4];

%Calcolo la Funzione obbiettivo: E(w)=(w1+w2*x(i)+w3*x(i)^2 -t(i))^2
Ew=0;
for i=1:n
    X=[1;xt(i);xt(i)^2;xt(i)^3];
    Ew=Ew+(W.'*X -t(i))^2;
end
Ew=Ew/2;
Erms=sqrt(2*Ew/n);

%Calcolo la derivata dEw/dw di Ew rispetto ai parametri w
DErms=jacobian(Erms,[w1 w2 w3 w4])

%Calcolo la soluzione analitica risolvendo dEw/dw=0
W_analitica=solve(DErms);
W1=[double(W_analitica.w1);double(W_analitica.w2);double(W_analitica.w3); double(W_analitica.w4);]

%Valuto la soluzione su tutti i punti x
for i=1:50
    X=[1;x(i);x(i)^2;x(i)^3];
    y_stimata(i)=W1.'*X;
end

plot(x,y_stimata,'g');