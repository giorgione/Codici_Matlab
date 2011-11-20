%Calcolo della Probabilità della Sequenza di Stati Visibili
%     t
% P( V ) attraverso l' algoritmo iterativo
%
%Considero la Catena di Markov a STATI NASCOSTI costituita da:
%
% - 4 STATI NASCOTI:  W0, W1, W2, W3
%
% - 4 STATI VISIBILI: V0,V1, V2, V3
%
% Matrice delle TRANSIZIONI DI STATO: A(i,j)= P( Wj(t)| Wi(t-1))
clc;clear;
A=[1   0    0   0;
   .2  .3  .1   .4;
   .2  .5  .2   .1;
   .8  .1  .0   .1];

% Matrice degli STATI VISIBILI: B(i,j)= P( Uj(t)| Ui(t) )
B=[1  0  0  0  0;
   0 .3 .4 .1 .2;
   0 .1 .1 .7 .1;
   0 .5 .2 .1 .2];

% Definisco la Sequenza di Stati Visibili ->5=T
V=[4 2 4 3 1];
%                                                                t 
%Calcolo tutte le possibili combinazioni degli Eventi: Tutte le C  possibili
% sequenze di STATI NASCOSTI W (osservabili dal modello HMM)--> cammini
% possibili di lunghezza 4
%
%genero tutte le combinazioni con ripetizione possibili degli indici degli
%stati : [1 2 3 4] --> 4^4=256
Vt=pick(1:4,4,'or');
[m,n]=size(Vt);

%
%              T                                               
%Calcolo  P( Wr ) la Probabilità della SEQUENZE di STATI NASCOSTI  di
%lunghezza r.
P_Wr_t=ones(m,1);


%            T    T                                               T       T
%Calcolo P( V | Wr ),Probabilità di osservare gli STATI VISIBILI V data Wr
%la SEQUENZA di STATI NASCOSTI  di lunghezza r.
P_Vt_Wr_t=ones(m,1);


P_Vt=0;
for i=1:m
    for j=1:n
        %calcolo la probabilità di osservare gli STATI NASCOSTI
        if j==1
            P_Wr_t(i)=1;
        else
            P_Wr_t(i)=P_Wr_t(i) * A(Vt(i,j),Vt(i,j-1));
        end
        
        %calcolo la probabilità di osservare gli STATI VISIBILI dati gli
        %Stati Nascosti
        P_Vt_Wr_t(i)=P_Vt_Wr_t(i)*B(Vt(i,j),V(j));
        
        %calcola la probabilità di osservare la sequenza di 
        P_Vt=P_Vt+P_Wr_t(i)*P_Vt_Wr_t(i);
    end
end

%Suppongo che lo stato iniziale e W1
Vt=pick(1:4,3,'or');
[m,n]=size(Vt);
Vt=[ones(m,1) Vt];
[m,n]=size(Vt);
%                                          T
%Calcolo la Probabilità della SEQUENZE P( Wr )
P_Wr_t=ones(m,1);
P_Vt_Wr_t=ones(m,1);
P_Vt=0;
for i=1:m
    for j=1:n
        if j==1
            J=2;
        else
            J=Vt(i,j-1);
        end
        I=Vt(i,j);
        K=V(j);
        %                                          T
        %Calcolo la Probabilità della SEQUENZE P( Wr )
        P_Wr_t(i)=P_Wr_t(i) * A(I,J);
        
        P_Vt_Wr_t(i)=P_Vt_Wr_t(i)*B(I,K);
        
        P_Vt=P_Vt+P_Wr_t(i)*P_Vt_Wr_t(i);
        
        display([ 'a(' num2str(I) ',' num2str(J) ')*b( ' num2str(I) ',' num2str(K) ')  ->' num2str(P_Wr_t(i)) ' * ' num2str(P_Vt_Wr_t(i)) ] )
       
        
    end
    V
    Vt(i,:)
end


