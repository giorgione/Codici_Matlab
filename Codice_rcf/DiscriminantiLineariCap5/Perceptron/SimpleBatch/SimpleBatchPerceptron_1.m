% SimpleBatchPerceptron_1.m
% 
% PROBLEMA DEL PERCETTRONE SEMLICE BIASED.
% Ricerca di un Discriminante lineare nello Spazio dei VETTORI AUGMENTED
% (vettore a) mediante la minimizazione della Funzione obbiettivo del 
% Percettrone:
%          
%           -----       t
% J(a) =    \       (-a) * y   -> a=[ao a1 a2]   ao=bias
%           /                     y=[1 x y]
%           -----
%           y appartente a Y --> Insieme dei Vettori missclassified
%
% Minimizzare questa funzione significa rendere minimo l'errore di
% Classificazione: scegliamo la configurazione dei parametri a tale che sia
% minimo il numero di pattern missclassified
%
% Procedimento:
%
% 1) Genero un insieme di Patter 2-D
%
% 2) Seleziono 10 Pattern casualmente e ne assegno 5 alla classe W1 , 5
%    alla classe W2
%
% 3) Effettuo il Labeling

clc;clear;close all;

k=menu('Seleziona i Dati:','1)Dati Non Separabili','2)Dati L.S. mediante retta non per l'' origine','3) Dati L.S. mediante retta per l''origine');
switch k
    case 1
        %Seleziona 10 vettori casualmente
        [x1 x2]=meshgrid(-10:.025:10);
        n=length(x1(:));
        indici=round(1+rand(1,10)*n);
        Y=[x1(:).';x2(:).'];
        Vet=Y(:,indici);  
        n=5;
    case 2
        % Il Metodo del Percettrone BIASED  converge alla Soluzione 
        % ottima poichè esiste una retta di Separazione passante per
        % l' origine che produce 0 errori di classificazione sul TRAINING
        % SET.
        x1=5+randn(1,10);
        y1=5+randn(1,10);

        x2=randn(1,10);
        y2=randn(1,10);
        
        %Insericsco i dati in un unico Vettore
        Vet=[x1 x2;y1 y2];
        n=length(x1);
        
    case 3
        %Il Metodo del Percettrone NOT BIASED converge alla soluzione OTTIMA
        % con 0 errori di classificazione sul TRAINING SET.
        x1=-2+randn(1,10);
        y1=-2+randn(1,10);

        x2=2+randn(1,10);
        y2=2+randn(1,10);
        
        %Insericsco i dati in un unico Vettore
        Vet=[x1 x2;y1 y2];
        n=length(x1);
end

%Suppongo che:
%- i primi 5 vettori (y1...y5) appartengono ad W1 --> label 1
%- gli altri 5 (y6..,y10) appartengano ad W2 --> label -1
V=[Vet(:,1:n) Vet(:,n+1:2*n)]
V=[ones(1,2*n);V]
V(:,n+1:2*n)=-V(:,n+1:2*n);
%Dati i 10 feature Vector voglio calcolare un discriminante lineare del tipo
%
% g(x)=a'.yi        tale che:
%
% a.'*yi >0   con yi=y1 e yi=y2
%
% a.'*yi <0   con yi=y3 e yi=y4
%
% Normalizzo il problema sostituendo tutti i campioni appartenti ad W2 con
% il loro negativi, semplificando il problema :
%Dati i 4 feature Vector voglio calcolare un discriminante lineare del tipo
%
% g(x)=a'.yi         tale che:
%
% a.'*yi >0   con yi=y1,y2,-y3,-y4
%
% Idea: Definisco una Funzione obbiettivo di a , J(a) che è minimizzata se
%       a è soluzione del problema..
%       Il calcolo della soluzione si riduce alla ricerca del minimo di una
%       funzione nello spazio dei PESI.
%       L' algoritmo di Ricerca è quello basato sul Gradiente.
%       La funzione obbiettivo va scelta oppotunamente

%Considero il criterio del Perceptron
%
% J(a)=sum(-a.'y) per ogni y in Ymissed
%
% Il gradiente di J(a) è:
%
% sum(-y) per ogni y in Ymissed
%
%Valore iniziale dei pesi a tale che: 
%       a.'*y <=0  -> errore di classificazione

options = optimset('Display','final','TolFun',1e-8,'MaxIter',500,'TolX',1e-8)
%}{

%Calcolo il minimo
[a ,Val,output]= fminsearch(@(x) Perceptron(x,V),randn(3,1),options)
plot3(a(1),a(2),a(3),'oy','MarkerSize',5,'MarkerFaceColor','y')


%Calcolo il prodotto e conto il numero di valori negativi che individua il
%numero di valori erroneamente classificati (Spazìo Augmented)
Errori=a.'*V;
NEr=length(find(Errori<0));



plot(Vet(1,1:n), Vet(2,1:n),'or','MarkerSize',5,'MarkerFaceColor','r');
hold on;
plot(Vet(1,n+1:2*n), Vet(2,n+1:2*n),'ob','MarkerSize',5,'MarkerFaceColor','b')


syms x y real;
S=solve(a.'*[1;x;y],y);
ezplot(S,[-10 10]);axis tight
title([ 'Retta di Separazione  numero errori:' num2str(NEr) ])

