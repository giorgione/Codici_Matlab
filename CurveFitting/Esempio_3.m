%Esempio di Curve Fitting Polinomiale
clc;clear;close all

x=linspace(0,1,150);
%curva da stimare
y=sin(2*pi*x);

%Valori osservati con Rumore additivo
y1=randperm(150);

%Numero punti osservati
n=20;
rumore=randn(1,n);
tt=y(y1(1:n))+rumore;
xt=x(y1(1:n));
figure(1);
plot(x,y,'r');hold on
plot(xt,tt,'ob');


%Approssimazione attraverso Polinomio di grado 3
basi=@(x)([1 x x^2 x^3 x^4 x^5 x^6 x^7 x^8 x^9]);
%calcolo i coefficienti del polinomio
[w,DesignMatrix]=LinearBasis(tt',xt,basi);

%Ipotizzo che  i dati osservati (ti,xi) siano caratterizzati da una distribuzione
%Gaussiana N(ti| y(xi,w), beta) con:
%
%  - media=y(xi,w)
%  - varianza= beta
%
% w rappresenta una variabile ALEATORIA di cui posso calcolare la
% distibuzione a posteriori mediante il teorema di bayes:
%               N
% Likelyhood= II   N(ti| y(xi,w), beta) --> funzione di (w,beta)
%               i=1
%
% Ipotizzo una prior su w --> P(w| a)=N(w| 0,alpha) funzione di (w,alpha)
% la media e la varianza sono noti
%
% Ottengo P(w|T,X,alpha,beta)
syms w0 w1 w2 w3 w4 w5 w6 w7 w8 w9

W=[w0;w1;w2; w3;w4; w5; w6; w7;  w8; w9];
M=length(W);

%HyperParametri di Precisione
Beta=11.1;
Alpha=0.0005;

Mo=zeros(M,1);
So=(1/Alpha)*eye(M);

%Utilizzo le formule dirette per calcolare La Posterior

SN_=inv(So)+Beta*DesignMatrix'*DesignMatrix;
SN=inv(SN_);
MN=SN*(inv(So)*Mo+Beta*DesignMatrix'*tt');

Posterior=GaussianaMulti(SN,MN,W);


%% Calcolo la stima dei parametri Massimizzando il log(Posterior)
DPosterior=jacobian(log(Posterior),[w0 w1 w2 w3 w4 w5 w6 w7 w8 w9])
W_analitica=solve(DPosterior);
W1=[double(W_analitica.w0);double(W_analitica.w1);double(W_analitica.w2);double(W_analitica.w3); double(W_analitica.w4);
    double(W_analitica.w5);double(W_analitica.w6);double(W_analitica.w7); double(W_analitica.w8) ; double(W_analitica.w9)];



%W0 è uguale a W1: i due approcci portano alla stessa soluzione

%% Calcolo w attraverso un problema di fitting con regolarizzazione
W=[w0; w1; w2; w3; w4; w5; w6; w7; w8 ;w9];
l=Alpha/Beta

%Calcolo la Funzione obbiettivo: E(w)=(w1+w2*x(i)+w3*x(i)^2 -t(i))^2
Ew=0;
for i=1:n
    X=basi(xt(i))';
    Ew=Ew+(W.'*X -tt(i))^2;
end
Ew=Ew/2  + l/2 *(W.'*W);

%Calcolo la derivata dEw/dw di Ew rispetto ai parametri w
DEw=jacobian(Ew,[w0 w1 w2 w3 w4 w5 w6 w7 w8 w9])

%Calcolo la soluzione analitica risolvendo dEw/dw=0
W_analitica=solve(DEw);
W2=[double(W_analitica.w0);double(W_analitica.w1);double(W_analitica.w2);double(W_analitica.w3); double(W_analitica.w4);
    double(W_analitica.w5);double(W_analitica.w6);double(W_analitica.w7); double(W_analitica.w8); double(W_analitica.w9)];


%W2 è uguale a W1  --> MASSIMIZZARE la LOG(Posterior)== Fitting SSD
%regolarizzato
%Valuto la soluzione su tutti i punti x
for i=1:150
    X=basi(x(i))'
    y_stimata(i)=W1.'*X;
end
figure(1)
plot(x,y_stimata,'m');hold on



%% Calcolo la Distribuzione Predittiva P(t |x,T,Alpha,Beta)
for i=1:150
  
    SigmaN=1/1+basi(x(i))*SN*basi(x(i))';
    MediaN=MN'*basi(x(i))';
    F(i)=  MediaN ;  %
    P(i)=GaussianaMulti(SigmaN,MediaN,x(i));
    Fup(i)= MediaN +SigmaN ; %GaussianaMulti(SigmaN,MediaN,x(i)+SigmaN);
    Fdown(i)=MediaN -SigmaN ; %GaussianaMulti(SigmaN,MediaN,x(i)-SigmaN); 
end
figure(1)
plot(x,F,'g')
plot(x,Fup,'y')
plot(x,Fdown,'y')

figure(2)
plot(x,P,'g')