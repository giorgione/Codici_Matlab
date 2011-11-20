%Calcola i Parent di 1 nodo X nel Modello Grafico
%
% parametri:
%
% X: nodo per il quale voglio ottenere i genitori
%
% Mxy: Matrice di aDiacenza dei Nodi tale che
%
%               |  1 se esiste il link     i ---> j P( xj | xi) 
%      M(i,j)= -|
%               |  0 se non esiste il link i ---->j P(xi | xj) 
%

function Parents=GetParents(X,Mxy)

%La colonna di Mxy(:,X) contiene tutti i possibili link verso X

Parents=find(Mxy(:,X)==1);
