% Descisione2Classi.m
% Problema di Decisione a 2 classi basato sul criterio MAP:
%
% C1: classe dei Salmoni
% 
% C2: classe dei Branzini
%
% La variabile aleatoria che utilizzo come caratteristica per la
% classificazione è la lungezza:
%
% Suppongo che le Verosimiglianze siano stimate mediante metodo di Parzen:
%
% P(x | C1)= N(1.5,0.9)  -> la classe C1 ha una Distribuzione della lunghezza
% di tipo Gaussiana con media 50 e Varianza 10
%
% P(x | C2)= N(4,.3)  -> la classe C1 ha una Distribuzione della lunghezza
% di tipo Gaussiana con media 50 e Varianza 10

clc;clear;close all;
%Estraggo i campioni
n=20;

xC1=3*rand(1,n);
xC2=rand(1,n);
plot(xC2,0,'ob','MarkerSize',5,'MarkerFaceColor','b');hold on;   
plot(xC1,0,'or','MarkerSize',5,'MarkerFaceColor','r');  drawnow 
hn=1/sqrt(n);
    
%Fisso il Volume 
Vn=hn;

%Applico il metodo di Parzen per la stima di P(x) con finestra Gaussiana
Px_c1=Parzen_Window_Simbolico(xC1,hn,Vn)

%Applico il metodo di Parzen per la stima di P(x) con finestra Gaussiana
Px_c2=Parzen_Window_Simbolico(xC2,hn,Vn)
figure;
h1=ezplot(Px_c1,[-5 5]);hold on;
set(h1,'LineWidth',2.2,'Color','r')

h2=ezplot(Px_c2,[-5 5]);
set(h2,'LineWidth',2.2,'Color','b');axis tight;drawnow

legend('P( X | C1 )','P( X | C2 )');
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
% Classifica C1 se P(x|c1)P(c1) > P(x|c2)P(c2)
% 
% Classifica C2 se P(x|c1)P(c1) < P(x|c2)P(c2)
%
% Basta quindi considerare la funzione Discriminante D(x):
%
%                 D(x) = P(x|c1)P(c1)-P(x|c2)P(c2)
%
% e studiarne la positività:
%
% D(x) > 0  classifica C1 -->Regione C1
%
% D(x) < 0  classifica C2 -->Regione C2
%

D= Px_c1*Pc1-Px_c2*Pc2;
figure;
h3=ezplot(D,[-5,10]);
set(h3,'LineWidth',2.2,'Color','b')
hold on;
title('Regioni di Decisione per Altezza: C1 se D(x)>0 C2 se D(x)<0')
legend('D(x): discriminante')
axis tight;drawnow