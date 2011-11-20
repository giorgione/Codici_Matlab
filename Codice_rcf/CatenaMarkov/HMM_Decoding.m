%function [P Alfa]=HMM_Decoding(A,B,V,StatoIniziale)
% Algoritmo HMM_Decoding (N^2*T) per calcolare la sequenza di stati nascosti
% aventi la massimma Probabilità di generare Sequenza di Stati
% Visibili V di lunghezza T:
%                              t
% P( ModelloMarkoviano(A,B) | V )
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
function [Path Alfa]=HMM_Decoding(A,B,V,StatoIniziale)

T=length(V);
[N,M]=size(A);
Alfa=zeros(T,N);

Path=zeros(1,T);
%Definisco le probabilità degli Stati Iniziali considerando il fatto ke lo
%stato iniziale sia W1 P(W=i)
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
   [Val J]=max(Alfa(t,:))
   Path(t)=J-1;
end

