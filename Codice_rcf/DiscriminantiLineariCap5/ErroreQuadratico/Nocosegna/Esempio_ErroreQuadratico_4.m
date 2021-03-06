% Esempio_ErroreQuadratico_1.m
% Ricerca di un Discriminante lineare nello Spazio dei Pesi (vettore a)
% Mi trovo in R2
clc;clear;close all;
k=menu('Seleziona i Dati:','1)Dati Non Separabili','2)Dati Linearmente Separabili');
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
        %DATI NON SEPARATI LINEARMENTE
        
        %estraggo i Dati per la Classe c1
        x1=EstrazioneNormale(20,2,10);
        y1=EstrazioneNormale(20,2,10);

        %estraggo i Dati per la Classe c2
        x2=EstrazioneNormale(1,2,10);
        y2=EstrazioneNormale(2,1,10);
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
%- i primi 2 vettori (x1,y1) appartengono ad W1 --> label 1
%- gli altri 2 (x2,y2) appartengano ad W2 --> label -1


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
% Idea: Definisco una Funzione obbiettivo di a , J(a) che � minimizzata se
%       a � soluzione del problema..
%       Il calcolo della soluzione si riduce alla ricerca del minimo di una
%       funzione nello spazio dei PESI.
%       L' algoritmo di Ricerca � quello basato sul Gradiente.
%       La funzione obbiettivo va scelta oppotunamente

%Considero il criterio del Perceptron
%
% J(a)=sum(-a.'y) per ogni y in Ymissed
%
% Il gradiente di J(a) �:
%
% sum(-y) per ogni y in Ymissed
%
%Valore iniziale dei pesi a tale che: 
%       a.'*y <=0  -> errore di classificazione

%[a1 a2]=meshgrid(-100:10:100);
%[m,n]=size(a1);
%Matrice m*n x 2
%a=[a1(:).'; a2(:).'];

%[J DJ]=ErroreQuadratico(a,V);
%J=reshape(J,m,n);
%surf(a1,a2,J)
%xlabel('a1');
%ylabel('a2');
%zlabel('J(a)');
%hold on;

%Calcolo il minimo
[a Val]= fminsearch(@(x) ErroreQuadratico(x,V),[1;1 ;-2],V)
plot3(a(1),a(2),Val,'oy','MarkerSize',5,'MarkerFaceColor','y')
ao=a(1);

[a1 a2]=meshgrid(-100:10:100);
[m,n]=size(a1);
%Matrice m*n x 2
a=[ao*ones(1,m*n);a1(:).'; a2(:).'];

[J DJ]=ErroreQuadratico(a,V);
J=reshape(J,m,n);
surf(a1,a2,J)
xlabel('a1');
ylabel('a2');
zlabel('J(a)');
hold on;



%learning rate fisso
lr=.05;
%tolleranza 
tol=0.0015;

ao=[ao 10 -21];
a=ao;
plot3(a(1),a(2),ErroreQuadratico(a.',V),'om','MarkerSize',5,'MarkerFaceColor','m')
k=1;
DJ=100;
while norm(lr*DJ,2)> tol && k<100
      
      %Calcolo il gradiente
      [J DJ]=ErroreQuadratico(a(k,:).',V);      
      lr*DJ
      %Applico il metodo Steepest Descendent
      a(k+1,:)=a(k,:)-lr*(DJ.');
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
%numero di valori erroneamente classificati (Spaz�o Augmented)
Errori=a.'*V;
NEr=length(find(Errori<0));



plot(Vet(1,1:n), Vet(2,1:n),'or','MarkerSize',5,'MarkerFaceColor','r');
hold on;
plot(Vet(1,n+1:2*n), Vet(2,n+1:2*n),'ob','MarkerSize',5,'MarkerFaceColor','b')


t=[-200;200];
syms x y;
S=solve(a(2:3).'*[x;y]+a(1))
ezplot(S,[-100 100])
title([ 'Retta di Separazione  numero errori:' num2str(NEr) ])
