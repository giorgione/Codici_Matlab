%Esempio di Curve Fitting Polinomiale
clc;clear;close all

x=linspace(0,1,50);
%curva da stimare
y=sin(2*pi*x);

%Valori osservati con Rumore additivo
y1=randperm(50);

%Numero punti osservati
n=10;
rumore=randn(1,n);
tt=y(y1(1:n))+rumore;
xt=x(y1(1:n));
figure(1);
plot(x,y,'r');hold on
plot(xt,tt,'ob');


%Approssimazione attraverso Polinomio di grado 3
basi=@(x)([1 x]);
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
syms w0 real 
syms w1 real 
%syms x  real 
%syms t real 
syms Alpha real 
syms Beta real
W=[w0;w1];
M=length(W);

Likelyhood=1;
for i=1:n
    Likelyhood=Likelyhood*GaussianaMulti(1/Beta,basi(xt(i))*W,tt(i));
end
L1=subs(Likelyhood,'Beta',1);
figure(2);
ezmesh(L1)

Mo=zeros(M,1);
So=(1/Alpha)*eye(M);
Prior=GaussianaMulti(So,Mo,W);
P1=subs(Prior,'Alpha',1);
figure(3);
ezmesh(P1)

Posterior=Likelyhood*Prior;
%Visualizzo la Posterior avendo Fissato gli IperParametri(Alpha,Betha):
%la Posterior è ancora una guassiana
figure(4);
Posterior1=L1*P1;
ezmesh(Posterior1)

%L'integrale rispetto a w non viene risolto: nn riesco a disegnare la
%Distribuzione Predittiva
%DistrPredittiva=int(int(GaussianaMulti(1,basi(x)*W,t)*Posterior1,w1),w0);

%Utilizzo le formule dirette

SN_=inv(So)+Beta*DesignMatrix'*DesignMatrix;
SN=inv(SN_);
MN=SN*(inv(So)*Mo+Beta*DesignMatrix'*tt');

Posterior2=GaussianaMulti(SN,MN,W);
P2=subs(Posterior2,{'Alpha','Beta'},{1,1});
figure(5);
ezmesh(P2)

%Calcolo la stima dei parametri Massimizzando il log(Posterior)
DPosterior=jacobian(log(Posterior1),[w0 w1])
W_analitica=solve(DPosterior);
W0=[double(W_analitica.w0);double(W_analitica.w1)]

%Calcolo la stima dei parametri Massimizzando il log(Posterior)
DPosterior=jacobian(log(P2),[w0 w1])
W_analitica=solve(DPosterior);
W1=[double(W_analitica.w0);double(W_analitica.w1)]

%W0 è uguale a W1: i due approcci portano alla stessa soluzione

%% Calcolo w attraverso un problema di fitting con regolarizzazione
W=[w0;w1;];
l=1;  %Alpha/Beta

%Calcolo la Funzione obbiettivo: E(w)=(w1+w2*x(i)+w3*x(i)^2 -t(i))^2
Ew=0;
for i=1:n
    X=[1;xt(i)];
    Ew=Ew+(W.'*X -tt(i))^2;
end
Ew=Ew/2  + l/2 *(W.'*W);

%Calcolo la derivata dEw/dw di Ew rispetto ai parametri w
DEw=jacobian(Ew,[w0 w1])

%Calcolo la soluzione analitica risolvendo dEw/dw=0
W_analitica=solve(DEw);
W2=[double(W_analitica.w0);double(W_analitica.w1)]

%W2 è uguale a W1  --> MASSIMIZZARE la LOG(Posterior)== Fitting SSD
%regolarizzato
%Valuto la soluzione su tutti i punti x
for i=1:50
    X=[1;x(i)];
    y_stimata(i)=W1.'*X;
end
figure(1)
plot(x,y_stimata,'m');hold on

MN=subs(MN,{'Alpha','Beta'},{1,1});
SN=subs(SN,{'Alpha','Beta'},{1,1});

%% Calcolo la Distribuzione Predittiva P(t |x,T,Alpha,Beta)
for i=1:50
  
    SigmaN=1/1+basi(x(i))*SN*basi(x(i))'
    MediaN=MN'*basi(x(i))'
    F(i)=GaussianaMulti(SigmaN,MediaN,x(i));
    Fup(i)=GaussianaMulti(SigmaN,MediaN,x(i)+SigmaN);
    Fdown(i)=GaussianaMulti(SigmaN,MediaN,x(i)-SigmaN); 
end
figure(1)
plot(x,F)
plot(x,Fup,'y')
plot(x,Fdown,'y')



