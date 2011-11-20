% Variable_Increment_Perceptron_Margin.m
%
% Ricerca di un Discriminante lineare nello Spazio dei Pesi (vettore a)
% mediante la procedura A CORREZIONE DI ERRORE CON MARGINE DEL SINGOLO PERCETTRONE AD 
% INCREMENTO VARIABILE
%          

% Procedimento:
%
% 1) Genero un insieme di Patter 2-D
%
% 2) Seleziono 10 Pattern casualmente e ne assegno 5 alla classe W1 , 5
%    alla classe W2
%
% 3) Effettuo il Labeling ed ottengo l' insieme V dei Pattern di Training
%
% 4) Applico la procedura A CORREZIONE DI ERRORE CON MARGINE DEL SINGOLO PERCETTRONE AD 
% INCREMENTO VARIABILE sull' insieme di TRAINING PATTERN.
% 
% 5) Disegno la retta di Separazione

clc;clear;close all;
n=10;

Vet=GeneraDati(2,n/2,'s')

%Suppongo che:
%- i primi 5 vettori (y1...y5) appartengono ad W1 --> label 1
%- gli altri 5 (y6..,y10) appartengano ad W2 --> label -1
V=[Vet(:,1:5) -Vet(:,6:10)];
plot(Vet(1,1:5), Vet(2,1:5),'or','MarkerSize',5,'MarkerFaceColor','r')
hold on;
plot(Vet(1,6:10), Vet(2,6:10),'ob','MarkerSize',5,'MarkerFaceColor','b')

%Mischio i Dati in maniera casuale
I=randperm(n);
Vet=Vet(:,I)

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% la procedura a correzione di Errore  con MARGINE del Singolo 
% Percettrone a incremento VARIABILE 


%Punto iniziale nello Spazio dei Parametri
ao=[1; -2];
a=ao;

plot(a(1),a(2),'om','MarkerSize',5,'MarkerFaceColor','m')


iter=1;
maxiter=500000;

%fisso il Margine
b=10;


%indice ultimo errore
prevbad = 1; 

%indice corrente
cidx = 2;

while cidx ~= prevbad && iter < maxiter
    
    %considero l'incremento n(k)= 1/k
    nk=1/cidx;
    
    %errore di classificazione per il pattern cidx
    if a.'*V(:,cidx)<= b
        %correggo il Vettore a
        a = a+nk*V(:,cidx);
        %aggiorna l' indice dell' ultimo pattern errato
        prevbad = cidx ;
    end
    %incrementa cidx: indice corrente
    cidx = mod(cidx , n) + 1;
    iter=iter+1;
end

%Calcolo il prodotto e conto il numero di valori negativi che individua il
%numero di valori erroneamente classificati
Errori=a.'*V;
NEr=length(find(Errori<0));


% t=-10:2:10;
% x=a(1)*t;
% y=a(2)*t;
% 
% plot(x,y,'g','LineWidth',2)
%
t=-100:50:500;
N=null(a.')
x=N(1)*t;
y=N(2)*t;

plot(x,y,'g','LineWidth',2)
title([ 'Simple-Sample Correction Retta di Separazione  numero errori:' num2str(NEr) ])


