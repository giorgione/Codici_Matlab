%function [P Beta]=HMM_Backward(A,B,V,StatoIniziale)
% Algoritmo HMM_Backward (N^2*T) per valutare la Probabilità della Sequenza di Stati
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
function [P Beta]=HMM_Backward(A,B,V,StatoIniziale)

T=length(V);
[N,M]=size(A);
Beta=zeros(T,N);

%Definisco le probabilità degli Stati Iniziali considerando il fatto ke lo
%stato iniziale sia W1 P(W=i)
P=zeros(1,M);
P(StatoIniziale+1)=1;

for j=1:N
    %Probabilità di osservare V(1) al tempo t=0 dato lo STATO WJ
    
    Beta(T,j)=1;%*B(j,k);  
end

for t=T-1:-1:1
   for i=1:N
        S=0;
        for j=1:N
            k=V(t+1)+1;
            S=S+A(i,j)*Beta(t+1,j)*B(j,k);
        end
        
        Beta(t,i)=S;
        
   end
end
Beta.'

P=Beta(1,StatoIniziale+1);
