%function [P Alfa]=HMM_Forward(A,B,V,StatoIniziale)
% Algoritmo HMM_Forward (N^2*T) per valutare la Probabilità della Sequenza di Stati
% Visibili V di lunghezza T:
%     t
% P( V  |  ModelloMarkoviano(A,B))
%
% dove:
%
%  V: sequenza di stati Visibili di lunghezza T
%
%  A: matrice N x N delle transizioni di stato (W0 W1 W2 ....Wn)
%
%  B: matrice degli stati T visibili (Vo V1 V2....Vt)
%
%  StatoIniziale: W 
function [P Alfa]=HMM_Forward(A,B,V,StatoIniziale)

T=length(V);
[N,M]=size(A);
Alfa=zeros(T,N);

%Definisco le probabilità degli Stati Iniziali considerando il fatto ke lo
%stato iniziale sia W(StatoIniziale+1)
P=zeros(1,M);
P(StatoIniziale+1)=1;
for j=1:N
    %Probabilità di osservare V(1) al tempo t=0 dato lo STATO WJ
    k=V(1)+1;
    Alfa(1,j)=P(j);%*B(j,k);  
end

for t=1:T-1
   for j=1:N
        S=0;
        for i=1:N
            S=S+Alfa(t,i)*A(i,j);
        end
        k=V(t+1)+1;
        Alfa(t+1,j)=S*B(j,k);
        
   end
end

P=0;
for i=1:N
    P=P+Alfa(T,i);
end