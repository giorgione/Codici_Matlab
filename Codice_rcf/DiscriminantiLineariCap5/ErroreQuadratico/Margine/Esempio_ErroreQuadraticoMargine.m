% Esempio_ErroreQuadraticoMargine.m
% Ricerca di un Discriminante lineare nello Spazio dei Pesi (vettore a)
% Mi trovo in R2
clc;clear;close all;

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


b=0.5;
%Calcolo il minimo
[a Val]= fminsearch(@(x) ErroreQuadraticoMargine(x,V,b),[1;1 ;-2],{V,b})
%plot3(a(2),a(3),Val,'oy','MarkerSize',5,'MarkerFaceColor','y')


%Calcolo il prodotto e conto il numero di valori negativi che individua il
%numero di valori erroneamente classificati
Errori=a.'*V;
NEr=length(find(Errori<0.5));

plot(Vet(1,1:n), Vet(2,1:n),'or','MarkerSize',5,'MarkerFaceColor','r');
hold on;
plot(Vet(1,n+1:2*n), Vet(2,n+1:2*n),'ob','MarkerSize',5,'MarkerFaceColor','b')


syms x y;
S=solve(a.'*[1;x;y],y);
ezplot(S,[-10 10]); axis tight;
title([ 'Retta di Separazione  numero errori:' num2str(NEr) ])