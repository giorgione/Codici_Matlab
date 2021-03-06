% Descisione2Classi.m
% Problema di Decisione a 2 classi basato sul criterio MAP:
%
% C1: classe dei Salmoni
% 
% C2: classe dei Branzini
%
% La variabile aleatoria che utilizzo come caratteristica per la
% classificazione � la lungezza:
%
% Suppongo che le Verosimiglianze siano stimate mediante metodo di Parzen:
%
% P(x | C1)= N(50,10)  -> la classe C1 ha una Distribuzione della lunghezza
% di tipo Gaussiana con media 50 e Varianza 10
clc;clear;close all;
%Estraggo i campioni
n=20;

xC1=5*rand(2,n);
xC2=3+2*rand(2,n);
    
hn=1/sqrt(n);
    
%Fisso il Volume 
Vn=hn;

%Applico il metodo di Parzen per la stima di P(x) con finestra Gaussiana
Px_c1=Parzen_2d(xC1,hn,Vn)

%Applico il metodo di Parzen per la stima di P(x) con finestra Gaussiana
Px_c2=Parzen_2d(xC2,hn,Vn)

ezsurf(Px_c1);
title('P( X | C1 )');

figure;
ezsurf(Px_c2);
title('P( X | C2 )');
%
% Supponiamo che la conoscenza a Priori sia:
%
% " Nel Mar Baltico il 70% della popolazione dei pesci � Salmone ed il 30%
% � Branzino "
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
% Classifica C1 se P(x|c1)P(c1) > P(x|c2)P(c2)
% 
% Classifica C2 se P(x|c1)P(c1) < P(x|c2)P(c2)
%
% Basta quindi considerare la funzione Discriminante D(x):
%
%                 D(x) = P(x|c1)P(c1)-P(x|c2)P(c2)
%
% e studiarne la positivit�:
%
% D(x) > 0  classifica C1 -->Regione C1
%
% D(x) < 0  classifica C2 -->Regione C2
%

D= Px_c1*Pc1-Px_c2*Pc2;
figure;
ezsurf(D);

hold on;
title('Regioni di Decisione per Altezza: C1 se D(x,y)>0 C2 se D(x,y)<0')
legend('D(x): discriminante')

figure;
ezcontour(D);hold on;
plot(xC1(1,:),xC1(2,:),'or')
plot(xC2(1,:),xC2(2,:),'ob')

