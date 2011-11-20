% Esempio_MSE_3.m
%
% Costruisco la Retta dei minimi quadrati per l' insieme di osservazioni
% applicando La risoluzione diretta del sistema di equazioni Normali:
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
% 6) Risolvo:
%                   t           t
%                  Y * Y * a = Y * b 
%
% 7) Disegno la Retta dei Minimi Quadrati


clc;clear;close all;
% carattere di newline
s=sprintf('\n');
n=10;

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
V=[Vet(:,1:n) Vet(:,n+1:2*n)];
V(:,n+1:2*n)=-V(:,n+1:2*n);


n=2*n;
%Mischio i Dati in maniera casuale
I=randperm(n);
V=V(:,I);


%Costruisco la Matrice delle Osservazioni
Y=V.';

%genero il vettore dei Margini b
b=.5*rand(n,1);

%Genero lo spazio dei Parametri a1,a2
[a1 a2]=meshgrid(linspace(-50,50,20));
[M,N]=size(a1);
a=[a1(:) a2(:)].';

figure;
%Calcolo la funzione del Percettrone e disegno
[J,DJ,H]=Mse(Y,a,b);
J=reshape(J,M,N);
surf(a1,a2,J)
xlabel('a1');
ylabel('a2');
zlabel('J(a)');
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algoritmo del gradiente per la ricerca del Minimo della funzione
% Percettrone

%Calcolo il minimo
[a Val]= fminsearch(@(x) Mse(Y,x,b),[1 ;-2],{Y,b});
plot3(a(1),a(2),Val,'oy','MarkerSize',5,'MarkerFaceColor','y')
display('Vettore dei Parametri della Retta Separatrice calcolati con fminsearch:')
disp(a)
disp(s)
Errori=a.'*V;
NEr=length(find(Errori<0));
display('Numero di Errori:');
disp(NEr)
disp(s)

%learning rate fisso
lr=.005;

%tolleranza
tol=1e-8;

%Punto iniziale nello Spazio dei Parametri
ao=[20 22];
a=ao;
plot3(a(1),a(2),Mse(Y,a.',b),'om','MarkerSize',5,'MarkerFaceColor','m')
k=1;
DJ=100;
while norm(lr*DJ,2)> tol && k<100
           
      %Calcolo il gradiente
      [J,DJ,H]=Mse(Y,a(k,:).',b); 
      
      lr=(DJ.'*DJ)/(DJ.'*H*DJ);
      % norm(lr*DJ,2)
      %Applico il metodo Steepest Descendent
      a(k+1,:)=a(k,:)-lr*(DJ.');
      k=k+1;
      
      %Valuto la funzione nel nuovo punto 
      [Jnew, DJnew, Hnew]=Mse(Y,a(k,:).',b);
      %Disgna il nuovo punto
      plot3(a(k,1),a(k,2),Jnew,'om','MarkerSize',5,'MarkerFaceColor','m')
      %Disegna la linea tra il nuovo ed il Vecchio Punto
      l=line( [a(k,1);a(k-1,1)] , [a(k,2);a(k-1,2)] , [Jnew ;J]);
      set(l,'Color','m','LineWidth',2);
      %pause
end

%Calcolo il prodotto e conto il numero di valori negativi che individua il
%numero di valori erroneamente classificati
a=a(end,:);
Errori=a*V;
NEr=length(find(Errori<0));
display('Vettore dei Parametri della Retta Separatrice calcolati con il mio gradiente:')
disp(a)
disp(s)
display('Numero di Errori:');
disp(NEr)
disp(s)

n=n/2;
figure;
plot(Vet(1,1:n), Vet(2,1:n),'or','MarkerSize',5,'MarkerFaceColor','r');
hold on;
plot(Vet(1,n+1:2*n), Vet(2,n+1:2*n),'ob','MarkerSize',5,'MarkerFaceColor','b')

%Calcola la linea di separazione   
syms x y;
S=solve(a*[x;y],y);
ezplot(S); axis tight;
title([ 'Retta di Separazione MSE  numero errori:' num2str(NEr) ])

