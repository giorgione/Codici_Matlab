% SingleSamplePerceptron.m
%
% Ricerca di un Discriminante lineare nello Spazio dei Pesi (vettore a)
% mediante la minimizazione della Funzione obbiettivo del Percettrone:
%          
%           -----
% J(a) =    \       (-a) * y
%           /
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
% 3) Effettuo il Labeling ed ottengo l' insieme V dei Pattern di Training
%
% 4) Applico l' Algoritmo del gradiente discendente per la Minimizazione
%    della funzione Percettrone per l' insieme di TRAINING PATTERN,
%    disegnando ad ogni passo lo spostamento lungo la superficie Soluzione

clc;clear;close all;
[x1 x2]=meshgrid(-10:.025:10);

%Seleziona 10 vettori 
%n=length(x1(:));
%indici=round(1+rand(1,10)*n);
%Y=[x1(:).';x2(:).'];
%Vet=Y(:,indici);

Vet=GeneraDati(2,5,'s')

%Suppongo che:
%- i primi 5 vettori (y1...y5) appartengono ad W1 --> label 1
%- gli altri 5 (y6..,y10) appartengano ad W2 --> label -1
V=[Vet(:,1:5) -Vet(:,6:10)];
plot(Vet(1,1:5), Vet(2,1:5),'or','MarkerSize',5,'MarkerFaceColor','r')
hold on;
plot(Vet(1,6:10), Vet(2,6:10),'ob','MarkerSize',5,'MarkerFaceColor','b')

%Mischio i Dati in maniera casuale
I=randperm(10);
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
% Algoritmo del gradiente per la ricerca del Minimo della funzione
% Percettrone


%Punto iniziale nello Spazio dei Parametri
ao=[1; -2];
a=ao;

plot(a(1),a(2),'om','MarkerSize',5,'MarkerFaceColor','m')

F=a.'*V;
[R C Val]=find(F<0)
k=-1;
maxiter=90000; 
iter=0;

while isempty(C)==0 && iter <=maxiter
    
       k=mod(k+1,10);
      
       %Il pattern k-esimo è Misclassfied
       if a.'*V(:,k+1) < 0
           %Aggiorna a
           a=a+V(:,k+1);
           
           %riclassifico i pattern
           %F=a.'*V;
           %[R C Val]=find(F<0);
       end
                 
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
t=-100:50:200;
N=null(a.')
x=N(1)*t;
y=N(2)*t;

plot(x,y,'g','LineWidth',2)
title([ 'Simple-Sample Correction Retta di Separazione  numero errori:' num2str(NEr) ])
