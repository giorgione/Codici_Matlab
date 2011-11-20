% RBF_Interpolazione_Regolarizzata.m
%Problema di Interpolazione mediante RBF REGOLARIZZATE: 
%                           mo
%Assegnati N punti x(i) in R   ed N numeri Reali d(i) voglio trovare la
%Funzione interpolante F tale che:
% 
%                            
% F(x(i))=d(i)   con i=1...N     Condizioni di INTERPOLAZIONE
%                                   MULTIDIMENSIONALE
%
% quindi:
% 
%           mo       1
%      F: R   ---> R    (INTERPOLAZIONE MULTIDIMENSIONALE)
%

%
% Gli RBF specificano il modello INTERPOLANTE:
%
%       ---
% F(x)= \   w(i) * teta(|| x - x(i) ||) = w(i)*teta(r) con r=||x -x(i)||
%       /
%       ----
%
% ovvero F è una combinazione lineare con pesi w(i) di Funzioni a Base Radiale 
% con centro x(i), generalmente sono funzioni non Lineari.
% Le Basi + utilizzate che mi consentono di ottenere un Sistema NON
% SINGOLARE sostituendo le condizioni di INTERPOLAZIONE sono:
%
% - Gaussiana:                     2      2
%                   teta(r)=exp(- r  / 2*o )      o>0
% 
% - Multiquadrica:                2    2   2
%                   teta(r)=sqrt(r  + o )/o
%  
%
% -Multiquarica Inversa:
%                            2        2    2  
%                   teta(r)=o / sqrt(r  + o )
%
% Cauchy:
%                            2    2    2  
%                   teta(r)=o / (r  + o )
% 
clc;clear;close all
k=menu('Seleziona il Modello Interpolante:','1)GAUSSIANA','2)MULTIQUADRICA','3)MULTIQUADRICA INVERSA','4)CAUCHY');
switch k
    case 1
        T=@(r,o) exp(-(r^2/o^2));
    case 2
        T=@(r,o)((r^2+o^2)^.5)/o^2;
    case 3
        T=@(r,o) o^2/((r^2+o^2)^.5);
    case 4
        T=@(r,o) o^2/(r^2+o^2);
end

%Problema di Interpolazione 2d su una Gaussiana

%Estraggo 10 punti da una Distribuzione Normale
X=randn(2,10);

Sigma=[1 .2;
       .2 1];

M=[0;0]
d=GaussianaMulti(Sigma,M,X).';
%Aggiungo Rumore a Casuale ai Dati;
Rumore=.2*randn(10,1);
display('Rumore')
disp(Rumore);
d=d+Rumore;

syms x1 x2;
F=GaussianaMulti(Sigma,M,[x1;x2]);
ezsurf(F);hold on;
title('Curva Originale e punti di interpolazione Rumorosi')
plot3(X(1,:),X(2,:),d+Rumore,'or','MarkerSize',4,'MarkerFaceColor','r')

%Costruico la Matrice di Interpolazione(matrice di Green)
o=1;
for i=1:10
    for j=1:10
        G(i,j)=T(norm(X(:,i)-X(:,j)),o);
    end
end


%calcolo i coefficienti dell' Interpolante
w=inv(G)*d;

%Genro i punti in cui andare a valutare la Funzione Interpolante
[x1 x2]=meshgrid(linspace(-3,3,40));
[m,n]=size(x1);
x=[x1(:).'; x2(:).'];
N=m*n;
for i=1:N
   %Valuro la Funzione interpolante nei punti
    Z1(i)=RBF_Interpolante(x(:,i),w,X,k,o);
    Z(i)=GaussianaMulti(Sigma,M,x(:,i)).';
end
display('Errore di Approssimazione:')
Errore=norm(Z-Z1);
disp(Errore)

Znew=reshape(Z1,m,n);
figure;
surf(x1,x2,Znew)
title('Superfice Interpolante non REGOLARIZZATA')


figure;
surf(x1,x2,(Znew-reshape(Z,m,n)).^2)
title(['Superfice Errore  errore globale:' num2str(Errore)])


lamda=input('Inserisci il Fattore di Regolarizzazione: ')
I=eye(10);

%calcolo i coefficienti dell' Interpolante
w=inv(G+lamda*I)*d;

%Genro i punti in cui andare a valutare la Funzione Interpolante
for i=1:N
   %Costruisco la Funzione interpolante
    Z2(i)=RBF_Interpolante(x(:,i),w,X,k,o);
end
display('Errore di Approssimazione:')
Errore=norm(Z-Z2);
disp(Errore)

Znew=reshape(Z2,m,n);
figure;
surf(x1,x2,Znew)
title('Superfice Interpolante REGOLARIZZATA')

figure;
surf(x1,x2,(Znew-reshape(Z,m,n)).^2)
title(['Superfice Errore  errore globale:' num2str(Errore)])

pause;
close all

