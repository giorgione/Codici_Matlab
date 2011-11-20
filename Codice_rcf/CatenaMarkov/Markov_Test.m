%Calcolo della Probabilità della Sequenza di Stati Visibili
%     t
% P( V ) attraverso :
%
% 1) l' algoritmo HMM Forward 
%
% 2) l' algoritmo HMM Backward
%
% 3) Risolvo il Problema di Decoding
%
%Considero la Catena di Markov costituita da:
%
% - 4 STATI NASCOTI:  W0, W1, W2, W3 
%
% - 1 STATO VISIBILE Iniziale assumibile  Vo 
%
% - 4 STATI VISIBILI: V1,V2, V3, V4
%
% Matrice delle TRANSIZIONI DI STATO: A(i,j)= P( Wj(t)| Wi(t-1))
clc;clear;

A=[   1   0   0    0;
     .2  .3  .1   .4;
     .2  .5  .2   .1;
     .8  .1  .0   .1];
%


% Matrice degli STATI VISIBILI: B(i,j)= P( Vj(t)| Wi(t) ) probabilità che
% lo stato NASCOSTO Wi assuma stato VISIBILE Vj
%  V0 V1 V2 V3 V4
B=[ 1  0  0  0  0 ;  %W0  
    0 .3 .4 .1 .2 ;  %W1
    0 .1 .1 .7 .1 ;  %W2
    0 .5 .2 .1 .2 ]; %W3 


% Definisco V la Sequenza di Stati Visibili ->5=T
%                      T
% Voglio calcolare p( V )
V=[3 1 3 2 0];

%Le Matrici sono indicizzate da 0

%Suppongo che lo stato iniziale sia W1 al tempo To=0
StatoIniziale=1;
P_StatoIniziale=1;

pause;
clc;
display('HHM Forward')
[P Alfa]=HMM_Forward(A,B,V,StatoIniziale)

pause;
clc;
display('HHM Backrward')
[P Beta]=HMM_Backward(A,B,V,StatoIniziale)

pause;
clc;
display('HHM Decoding')
[Path Alfa]=HMM_Decoding(A,B,V,StatoIniziale)