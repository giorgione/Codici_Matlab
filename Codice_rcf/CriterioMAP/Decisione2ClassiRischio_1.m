% Problema di Decisione a 2 classi basato sul criterio MAP,
% è nota anche la matrice dei Costi A(2 x 2):
%
% A(i,j): costo associato nell' effettuare l'azione a(i) quando lo stato in
%         natura è C(j).
%
% Suppongo che ci siano solo 2 azioni:
%
%  a1 -> decido che lo stato in natura è C1 
%  a2 -> decido che lo stato in natura è C2
%
% Definisco la matrice 0-1 di perdita, essa associa:
%
% - costo 0 ad una decisione correta    ->     A(1,1)=0=A(2,2)
%
% - costo 1 ad una decisione sbagliata  ->     A(1,2)=1=A(2,1)
%
clear;clc;close all
A=[0 .5;
   .25 0];
%
%
% C1: classe dei Salmoni
% 
% C2: classe dei Branzini
%
% La variabile aleatoria che utilizzo come caratteristica per la
% classificazione è la lungezza:
%
% Suppongo che le Verosimiglianze siano del tipo:
%
% P(x | C1)= N(50,10)  -> la classe C1 ha una Distribuzione della lunghezza
%                         di tipo Gaussiana con media 50 e Varianza 10

syms x ;
X=[x];
M=50;
%Matrice di Varianza Covarianza
S=[10^2 ];
Px_c1=GaussianaMulti(S,M,X);

%
% P(x | C2)= N(80,5)  -> la classe C1 ha una Distribuzione della lunghezza
%                          di tipo Gaussiana con media 80 e Varianza 5
M=[60];
%Matrice di Varianza Covarianza
S=[5^2 ];
Px_c2=GaussianaMulti(S,M,X);

figure;
h1=ezplot(Px_c1,[0 100]);hold on;
set(h1,'LineWidth',2.2,'Color','r')
h2=ezplot(Px_c2,[0 100]);
set(h2,'LineWidth',2.2,'Color','b')
%
% Supponiamo che la conoscenza a Priori sia:
%
% " Nel Mar Baltico il 70% della popolazione dei pesci è Salmone ed il 30%
% è Branzino "
%
% che si traduce in:
%
% P(C1)=7/10
%
% P(C2)=3/10
Pc1=0.7;
Pc2=0.3;

% Il criterio MAP classifica nel seguente modo:
%  
% Classifica C1 se (A(2,1)-A(1,1))*P(x|c1)P(c1) > (A(1,2)-A(2,2))*P(x|c2)P(c2)
% 
% Classifica C2 se (A(2,1)-A(1,1))*P(x|c1)P(c1) < (A(1,2)-A(2,2))*P(x|c2)P(c2)
%
% Basta quindi considerare la funzione Discriminante D(x):
%
%                 D(x) = (A(2,1)-A(1,1))*P(x|c1)P(c1) - (A(1,2)-A(2,2))*P(x|c2)P(c2)
%
% e studiarne la positività:
%
% D(x) > 0  classifica C1 -->Regione C1
%
% D(x) < 0  classifica C2 -->Regione C2
%
%
% Nel caso della Matrice A 1-0 loss ottengo la regola di decisione Baysiana
% semplice essendo le differenze dei termini di A uguali ad 1;
%
%       D(x) = P(x|c1)P(c1)-P(x|c2)P(c2)
%
% per cui:
%
% Minimizzare il RISCHIO equivale alla regola di decisione Bayesiana

D= (A(2,1)-A(1,1))*Px_c1*Pc1 - (A(1,2)-A(2,2))*Px_c2*Pc2;
figure;
h3=ezplot(D,[0,100]);
set(h3,'LineWidth',2.2,'Color','b')
hold on;
h4=line([0 100],[0 0]);
set(h4,'LineWidth',1.2,'Color','r')
title('Regioni di Decisione per Altezza: C1 se D(x)>0 C2 se D(x)<0')

Zero=solve(D);
plot(double(Zero),[0;0],'bo','MarkerFaceColor','r');

%Calcolo la Probabilità di Errore (assume il minimo Valore Possibile):

%integrale di P(x|c2)P(c2) nella Regione R1
P2=double (int((A(1,2)-A(2,2))*Px_c2*Pc2,x,-inf,Zero(1))+int((A(1,2)-A(2,2))*Px_c2*Pc2,x,Zero(2),inf));

%integrale di P(x|c1)P(c1) nella Regione R2
P1=double (int((A(2,1)-A(1,1))*Px_c1*Pc1,x,Zero(1),Zero(2)));