% Esempio_ErroreQuadratico.m
%
% Ricerca di un Discriminante lineare nello Spazio dei Pesi (vettore a)
% Mi trovo in R2
% PROBLEMA DEL PERCETTRONE SEMLICE NOT BIASED.
% Ricerca di un Discriminante lineare nello Spazio dei Pesi (vettore a)
% mediante la minimizazione della Funzione obbiettivo dell' Errore
% Quadratico: 
%          
%           -----        t     2
% J(a) =    \      ( (a) * y )     -> a=[a1 a2]
%           /                          y=[x y]
%           -----
%           y appartente a Y --> Insieme dei Vettori missclassified
%
% NB: Ricerco la retta di Separazione migliore passante per l' origine
% attraverso l' ALGORITMO DEL GRADIENTE DISCENDENTE IMPLEMENTATO DA ME...
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
% carattere di newline
s=sprintf('\n');

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
V=[Vet(:,1:n) -Vet(:,n+1:2*n)]


%Genero lo spazio dei Parametri
[a1 a2]=meshgrid(linspace(-100,100,20));
[M,N]=size(a1);
a=[a1(:) a2(:)];


%Calcolo la funzione del Percettrone e disegno
[J DJ]=ErroreQuadratico(a.',V);
J=reshape(J,M,N);
surf(a1,a2,J)
xlabel('a1');
ylabel('a2');
zlabel('J(a)');
hold on;



%Dati i 4 feature Vector voglio calcolare un discriminante lineare del tipo
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



%Calcolo il minimo
[a Val]= fminsearch(@(x) ErroreQuadratico(x,V),[1;1],V);
plot3(a(1),a(2),Val,'oy','MarkerSize',5,'MarkerFaceColor','y')
display('Vettore dei Parametri della Retta Separatrice calcolati con fminsearch:')
disp(a)
disp(s)
Errori=a.'*V;
NEr=length(find(Errori<0));
display('Numero di Errori:');
disp(NEr)
disp(s)

%Genero lo spazio dei Parametri
[a1 a2]=meshgrid(linspace(-5000,5000,20));
[M,N]=size(a1);
a=[a1(:) a2(:)];


%Calcolo la funzione del Percettrone e disegno
[J DJ]=ErroreQuadratico(a.',V);
J=reshape(J,M,N);
surf(a1,a2,J)
xlabel('a1');
ylabel('a2');
zlabel('J(a)');
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algoritmo del gradiente per la ricerca del Minimo della funzione
% Percettrone

%learning rate fisso
lr=0.05;

%tolleranza
tol=1e-10;

%Punto iniziale nello Spazio dei Parametri
ao=[300 200];
a=ao;
plot3(a(1),a(2),ErroreQuadratico(a.',V),'om','MarkerSize',5,'MarkerFaceColor','m')
k=1;
DJ=100;
while norm(lr*DJ,2)> tol && k<100
      
      %Calcolo il gradiente
      [J DJ]=ErroreQuadratico(a(k,:).',V);
      
      %Applico il metodo Steepest Descendent
      a(k+1,:)=a(k,:)-lr*(DJ.');
      k=k+1;
      %lr=1/k;
      
      %Valuto la funzione nel nuovo punto 
      [Jnew DJnew]=ErroreQuadratico(a(k,:).',V);
      %Disgna il nuovo punto
      plot3(a(k,1),a(k,2),Jnew,'om','MarkerSize',5,'MarkerFaceColor','m')
      %Disegna la linea tra il nuovo ed il Vecchio Punto
      l=line( [a(k,1);a(k-1,1)] , [a(k,2);a(k-1,2)] , [Jnew ;J]);
      set(l,'Color','m','LineWidth',2);
      pause
      
      
      
end

A=a(end,:);


%Calcolo il prodotto e conto il numero di valori negativi che individua il
%numero di valori erroneamente classificati (Spazìo Augmented)
Errori=A(end,:)*V;
NEr=length(find(Errori<0));
display('Vettore dei Parametri della Retta Separatrice calcolati con il Gradiente Discendente:')
disp(A)
disp(s)
NEr=length(find(Errori<0));
display('Numero di Errori:');
disp(NEr)
disp(s)

figure;
plot(Vet(1,1:n), Vet(2,1:n),'or','MarkerSize',5,'MarkerFaceColor','r');
hold on;
plot(Vet(1,n+1:2*n), Vet(2,n+1:2*n),'ob','MarkerSize',5,'MarkerFaceColor','b')


t=[-10;10];
syms x y;
S=solve(A*[x;y],y);axis tight;
ezplot(S,[-10 20])
title([ 'Retta di Separazione  numero errori:' num2str(NEr) ])
