%Esercizio 1
%
%Tecnica di campionamento INVERSO (CONTINUO) sulla Distribuzione esponenziale
%
%                           -x/lambda
%Pdf(x| lamda)= 1/lamda *  e
%                 
%
%                   -x/lambda
%Cdf(x| lamda)=1 - e
%
%Nb: La distribuzione esponenziale è una Distribuzione di Probabilità
%   CONTIUNA che descrive la DURATA DI VITA di un fenomeno che non invecchia
%   - PRIVO DI MEMORIA. Es: decadimento di una particella RadioAttiva
%
clc;close all;clear
PDFesp=@(x,lambda) (1/lambda)*exp(-x./lambda);
CDFesp=@(x,lambda) 1 - exp(-x./lambda);
InvCDFesp=@(U,lambda) -log(1-U)*lambda;
K = 1000;

X=linspace(0,10,K);
lambda=0.5;

Y=PDFesp(X,lambda);
Y1 = exppdf(X,lambda);
 
 % create a new figure
figure( 1 ); clf;
subplot(2,1,1)
plot(X,Y);
title('My exppdf')
subplot(2,1,2)
plot(X,Y1) ;
title('exppdf')

%Campiono la PDF con le builtin functions
Y = exprnd(lambda,1,K);
bins=1:10;

% Show the histogram of the simulated draws
counts = hist( Y , bins );
figure( 2 ); clf;
subplot(2,1,1)
bar(  bins , counts , 'k' );
xlim( [ 0 10 ] );
xlabel( 'x' );
ylabel( '(1/lambda)*exp(-x./lambda) ' );
title( 'Campionamento Simulato sulla distr. Esponenziale' );




%genero K campioni Random a distribuzione uniforme in [0 1]
Y1=unifrnd(0,1,1,K);
outcome=zeros(1,K);
%Campionamento Inverso Continuo

outcome=InvCDFesp(Y1,lambda);


% Show the histogram of the simulated draws
counts = hist( outcome , bins );
subplot(2,1,2)
bar( bins , counts , 'k' );
xlim( [ 0  10 ] );
xlabel( 'x' );
ylabel( '(1/lambda)*exp(-x./lambda) ' );
title( 'Campionamento Simulato sulla distr. Esponenziale' );