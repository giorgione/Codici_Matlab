% Esempio_MSE_2.m
%
% Costruisco la Retta dei minimi quadrati per l' insieme di osservazioni
% applicando l' algoritmo del Gradiente alla funzione MSE:
%
%   t           t
%  Y * Y * a = Y * b 
%
%
% Procedimento:
%
% 1) Genero un Set di Dati 2-D
%
% 2) Seleziono 10 Pattern casualmente e ne assegno 5 alla classe W1 , 5
%    alla classe W2
%
% 3) Effettuo il Labeling ed ottengo l' insieme V dei Pattern di Training
%
% 4) Costruisco la Matrice delle osservazioni Y (10 x 2)
%
% 5) Costruisco il Vettore dei Margini b (10 x 1)
%
% 6) Applico il Metodo del Gradiente per il calcolo del Vettore dei
%    Coefficienti
%
% 7) Disegno la Retta dei Minimi Quadrati


clc;clear;close all;
n=10;

k=menu('Seleziona i Dati:','1)Dati Non Separabili','2)Dati L.S. mediante retta non per l'' origine','3) Dati L.S. mediante retta per l''origine');
switch k
    case 1
        %DATI NON SEPARATI LINEARMENTE
        
        %Seleziona 10 vettori casualmente
        [x1 x2]=meshgrid(-10:.025:10);
        n=length(x1(:));
        indici=round(1+rand(1,10)*n);
        Y=[x1(:).';x2(:).'];
        Vet=Y(:,indici);  
        n=5;
     case 2
        % Il Metodo del Percettrone NOT BIASED non converge alla Soluzione 
        % ottima poichè non esiste una retta di Separazione passante per
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

%Suppongo che:
%- i primi 5 vettori (y1...y5) appartengono ad W1 --> label 1
%- gli altri 5 (y6..,y10) appartengano ad W2 --> label -1
plot(Vet(1,1:n), Vet(2,1:n),'or','MarkerSize',5,'MarkerFaceColor','r');
hold on;
plot(Vet(1,n+1:2*n), Vet(2,n+1:2*n),'ob','MarkerSize',5,'MarkerFaceColor','b')




n=2*n;
%Mischio i Dati in maniera casuale
I=randperm(n);
Vet=Vet(:,I)

%Costruisco la Matrice delle Osservazioni
Y=V.';

%genero il vettore dei Margini b
b=.5*rand(n,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algoritmo del gradiente per la ricerca del Minimo della funzione
% Percettrone

%Calcolo il minimo
[a Val]= fminsearch(@(x) Mse(Y,x,b),[0;1 ;-2],{Y,b})
plot3(a(2),a(3),Val,'oy','MarkerSize',5,'MarkerFaceColor','y')


%Calcolo il prodotto e conto il numero di valori negativi che individua il
%numero di valori erroneamente classificati
Errori=a.'*V;
NEr=length(find(Errori<0));



t=[-200;200];
syms x y;
S=solve(a(2:3).'*[x;y]+a(1),y);
ezplot(S,[-100 100]); axis tight;
title([ 'Retta di Separazione MSE  numero errori:' num2str(NEr) ])
