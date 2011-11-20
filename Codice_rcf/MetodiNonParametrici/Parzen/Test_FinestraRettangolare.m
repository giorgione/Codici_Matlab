% Test sulla funzione Rettangolo:
%
% Voglio Verificare che La Stima di Parzen di una distribuzione di
% Probabilit� ha integrale uguale ad 1

clear;clc;close all;
n=10;

x=rand(1,n);


%Rettangolo centrato in 5 di base (Ampiezza) 2*12.5=25 e Altezza 1: Area=25
ezplot('rettangolo((x-5)/5,5)',[-20 20])
Vn=quad(@(x) rettangolo((x-5)/5,5),-20,20)

%Rettangolo centrato in 5 di base (Ampiezza) 2*2.5=5 e Altezza 1: Area=5
figure;
ezplot('rettangolo((x-5)/5,1)',[-20 20])
Vn1=quad(@(x) rettangolo((x-5)/5,1),-20,20)

%Verifico che l'integrale di rettangolo((x-5)/5,5)/Vn � uguale ad 1:
% ho normalizzato
quad(@(x) rettangolo((x-5)/5,5)/Vn,-20,20)

quad(@(x) rettangolo((x-5)/5,1)/5,-20,20)

hn=5;
%Verifico che l'integrale di Pn � 1: Pn � una distribuzione di Probabilit�
Pn=@(X) SommaRettangoli2(X,x,hn)/n
I=quad(Pn,-20,20)

figure;
ezplot(Pn,[-20 20])