%Stima della Distribuzione Predittiva a Priori sul Modello grafico del
%Consenso caratterizzato dalle seguenti Variabili Aleatorie
%
% u:  Variabile che rappresenta la risposta ad una domanda del tipo
%     "quanti abitanti ci sono ad Irvine California" fornita ad un
%     insieme di M Individui per N volte. Ogni volta un indivuo
%     fornisce una risposta y(i,j).
%
% Problema: stimare il valore esatto di u
%
% Variabili presenti nel MODELLO GRAFICO:
%
% u: sconosciuto - latente
%
% y(i,j) ~ Norm( y(i,j)|u,v(j) ) --> Ogni utente da una stima di u con una 
%                                    certa precisione v(j)
% i=1...M  --> numero di Individui
% j=1...N  --> numero di Stime effettuat per ciauscun individuo
%
% v(j) = 1/tau(j)  --> variabile DETERMINISTICA
%
% tau(j) ~ Gamma( tau(j) |a,b) --> Precisione della Gaussiana modellata con
%                                   Gamma distr. di param (a,b)
%
clc;clear;close all;
% Supponiamo che il valore di u sia 150
% Stimare la distribuzione predittiva p(y(i,j))
% Dati Noti
M=3; %Individui
N=10; %Risposte fornite nel tempo

% (a,b,u) sono parametri noti ma potrebbere essere sostituiti con delle prior
u=150;
a=0.3; 
b=0.33; 
%Numero di simulazioni
T=1;
%Dati Sconosciuti
% Y(i,j) = Risposte fornite dagli utenti
%
% Tau(kj) = Precisione degli Utenti

%
% La Distribuzione congiunta del Modello Grafico è 
%
% P(X) = P(y(i,j),tau(j),)
%
%                    M                     N
%                   T T                   T T
%   p(a)*p(b)*p(u)  | | Gamma(tau(j),a,b) | | Norm(y(i,j)|u,v(j))
%   --------------  j=1                   i=1
%        ^
%        |_ vanno via perchè noti
%
% 
%Inizializzo tutti le Variabili presenti nel Modello
NVar=M+M*N;
P=@(Theta)WrapperJointConsensus(Theta,M,N,u,a,b);

%INIZIALIZZO I CAMPIONI 
%Campioni che genero nel tempo: ogni colonna rappresenta il campione
%estratto dal Processo di campionamento
Theta=zeros(NVar,T);
seed=1; rand( 'state' , seed ); randn('state',seed ); % set the random seed

Theta(1:3,1)=3*rand(3,1); %Varianza Iniziale degli Individui
Theta(4:end,1)=u+3*rand(M*N,1); %Scelte y(i,j)

%Campiono
Samples=Fun_MetroPolisHastingsSampler_CW(P,Theta,NVar,T,[2]);

yij=Samples(4:end,:);

%Ogni colonna di yij rappresenta j(i,j) al passo t
Individui=1:M;
Individui=repmat(Individui,1,N);
Individui=repmat(Individui.',1,T);
plot(Individui.',yij.','or')

%Calcolo p(u,tau|yij,a,b)
%
% yij --> campioni generati in precedenza
NVar=M+1;
T=500;
yij=reshape(yij,M,N);

P=@(Theta)WrapperJointConsensus2(Theta,M,N,yij,a,b);
Theta=zeros(NVar,T);
Theta(1:3,1)=rand(3,1); %Varianza Iniziale degli Individui
Theta(4:end,1)=u+10*rand(1,1); %Scelte u
Samples=Fun_MetroPolisHastingsSampler_CW(P,Theta,NVar,T,[2]);
figure(2)
plot(1:T,Samples(4,:),'ro-')
xlabel('u')
ylabel('p(u)')
figure(3)
plot(1:T,Samples(1,:),'ro-');hold on
plot(1:T,Samples(2,:),'bo-');hold on
plot(1:T,Samples(3,:),'go-');hold on
xlabel('thau')
ylabel('p(thau)')

