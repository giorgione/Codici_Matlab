% EM_Stima_1.m
%
% Stima dei parametri di un Modello Gaussiano attraverso l' Algoritmo EM
% supponendo che i dati siano Mancanti
%
clc;clear;close all;
%siamo in R2
Teta0=[0;0;1;1];
syms x4 real

% Matrice delle osservazioni (i-esima colonna --> osservazione xi )
Data=[0 1 2 x4;
      2 0 2 4];
D=double(Data(:,1:3));
%Disegno i 3 punti
plot(D(1,1),D(2,1),'or','MarkerFaceColor','r');hold on;
plot(D(1,2),D(2,2),'ob','MarkerFaceColor','b');
plot(D(1,3),D(2,3),'og','MarkerFaceColor','g');


%Costruisco la funzione di Verosimiglianza dei dati completi
Q=@(t) (-1)* ( GaussianaMulti_Punti(diag(t(3:4).^2),[t(1);t(2)],double(Data(:,1:3)))   - (1+t(1).^2)/(2*t(3).^2) -((4-t(2))^2)/(2*t(4)^2)-log(2*pi*t(3)*t(4)));

%M-Step
options = optimset('Display','final','MaxFunEvals',500,'MaxIter',500,'TolFun',1e-6,'TolX',1e-6);
T(:,1)=Teta0;
k=1;
while (k<8)
    
    %Disegna l' Ellisse definita dai parametri stimati
    [x,y]=DisegnaEllisse(T(1:2,k),T(3,k),T(4,k),40,0);
    plot(x,y,'Color',[1/k 1 0]);hold on
    
    [T(:,k+1),Val,exitflag,output] = fminsearch(Q,[T(1:2,k);sqrt(T(3:4,k))],options);
    T(3:4,k+1)=T(3:4,k+1).^2;
    
    pause;
    k=k+1;
end

