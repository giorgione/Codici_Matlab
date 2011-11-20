%
% Ricerca di un Discriminante lineare nello Spazio dei Pesi (vettore a)
% Mi trovo in R2
clc;clear;close all;
[x1 x2]=meshgrid(-10:0.25:10);

%Seleziona 10 vettori 
n=length(x1(:));
indici=round(1+rand(1,10)*n);
Y=[x1(:).';x2(:).'];
Vet=Y(:,indici);

V=[Vet(:,1:5) -Vet(:,6:10)]
%Suppongo che:
%- i primi 2 vettori (y1,y2) appartengono ad W1 --> label 1
%- gli altri 2 (y3,y4) appartengano ad W2 --> label -1


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

[a1 a2]=meshgrid(-100:10:100);
[m,n]=size(a1);
%Matrice m*n x 2
a=[a1(:).'; a2(:).'];

[J DJ]=ErroreQuadratico(a,V);
J=reshape(J,m,n);
surf(a1,a2,J)
xlabel('a1');
ylabel('a2');
zlabel('J(a)');
hold on;

%Calcolo il minimo
[a Val]= fminsearch(@(x) ErroreQuadratico(x,V),[1 ;-2],V)
plot3(a(1),a(2),Val,'oy','MarkerSize',5,'MarkerFaceColor','y')

%learning rate fisso
lr=.05;
%tolleranza
tol=0.015;

ao=[1 -2];
a=ao;
plot3(a(1),a(2),ErroreQuadratico(a.',V),'om','MarkerSize',5,'MarkerFaceColor','m')
k=1;
DJ=100;
while norm(lr*DJ,2)> tol && k<100
      
      %Calcolo il gradiente
      [J DJ]=ErroreQuadratico(a(k,:).',V);      
      % norm(lr*DJ,2)
      %Applico il metodo Steepest Descendent
      a(k+1,:)=a(k,:)-lr*(DJ.')
      k=k+1;
      
      %Valuto la funzione nel nuovo punto 
      [Jnew DJnew]=ErroreQuadratico(a(k,:).',V);
      %Disgna il nuovo punto
      plot3(a(k,1),a(k,2),Jnew,'om','MarkerSize',5,'MarkerFaceColor','m')
      %Disegna la linea tra il nuovo ed il Vecchio Punto
      l=line( [a(k,1);a(k-1,1)] , [a(k,2);a(k-1,2)] , [Jnew ;J]);
      set(l,'Color','m','LineWidth',2);
      pause
end

%Calcolo il prodotto e conto il numero di valori negativi che individua il
%numero di valori erroneamente classificati
Errori=a.'*V;
NEr=length(find(Errori<0));

%Calcola la linea di separazione   
N=null(a.');
m=N(1)/N(2);
t=-10:10;
y=m*t;
figure;
plot(t,y,'g','LineWidth',2)
title([ 'Funzione Errore Quadratico - Retta di Separazione  numero errori:' num2str(NEr) ])
hold on

plot(Vet(1,1:5), Vet(2,1:5),'or','MarkerSize',5,'MarkerFaceColor','r')
hold on;
plot(Vet(1,6:10), Vet(2,6:10),'ob','MarkerSize',5,'MarkerFaceColor','b')