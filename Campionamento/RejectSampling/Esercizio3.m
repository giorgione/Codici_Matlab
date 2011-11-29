% Esercizio 3 
%
% CAMPIONAMENTO CON RIGETTO su una distribuzione esponenziale
%
%                           -x/lambda
% Pdf(x| lamda)= 1/lamda *  e              --> TARGET DISTR.

clc;clear;close all
%Considero N punti
N=100;
syms x

%Definizione trovata su Bishop
P=BetaDistribution(x,2,1);

%Definizione trovata sul testo degli esercizi
Beta=@(x,a,b)((x.^(a-1)).*(1-x).^(b-1))./beta(a,b);

%Probabilità Beta(x,2,1) è una pdf per 1<x<0
P=Beta(x,2,1)
%Verifico che è una distribuzine per 1<x<0
int(P,0,1);
X=linspace(0,1,N);
plot(X,Beta(X,2,1)); hold on;

%Considero un UpperBound per P(x) qualsiasi
%
% P1= 2*P(x)

P11=2*P;
P1=P11/int(P11,0,1); %Normalizzo
int(P1,0,1)


PDFUnif=ones(1,N)/N;
plot(X,PDFUnif,'r');

%Calcolo il coefficiente C tale che
%PDFUnif(x) * c > P(x) per ogni x
maxP=Beta(1,2,1);
c=maxP*N;

plot(X,c*PDFUnif,'g');

CDFUnif=cumsum(PDFUnif)
%Campiono mediante campionamento inverso
K=100;
outcome=zeros(1,K);
Y1=unifrnd(0,1,1,K);

%Campionamento Inverso
for i=1:K
    j=1;
    while Y1(i) >= CDFUnif(j)
        j=j+1;          
    end
    outcome(i)=(j-1)/N;
end

%genero K campioni Random dalla distribuzione uniforme in [0 c/N]

Y1=unifrnd(0,c/N,1,K);

j=1;
for i=1:K
    
    if(Y1(i) < Beta(outcome(i),2,1))
       campioni(j)=outcome(i); 
       plot(campioni(j),Y1(i),'om')
       j=j+1;
    end
   
end