%Test distribuzione Circolare Uniforme
%           1/pi*r^2    x^2 + y^2 < R^2
% P(x,y) =  0
%
% 
% va convertita in forma polare
%
%           r/pi*r^2     0 <  r < Raggio  U 0 <  t < 2*pi
% P(r,t) =  0           

%Calcolo l'area di una circonferenza
clc;clear;close all
syms r
F=r;
Raggio=4;
A=int(int(F,'r',0,Raggio),'t',0,2*pi)

%Verifico che la pdf sommi ad 1 
A=int(int(r/A,'r',0,Raggio),'t',0,2*pi)
 
%Cumulativa rispetto gli estremi di Integrazione r:[0 R] t:[0 R]
syms R a
A=int(int(r/A,'r',0,R),'t',0,a)

%Campionare una distribuzione del genere e TROPPO COMPLICATO
T=1500;

%Utilizzo MCMC
%Inizializzo tutti le Variabili presenti nel Modello
NVar=2;
P=@(Theta)DistrCircolare(Theta,[10;10],4);
 

%INIZIALIZZO I CAMPIONI 
%Campioni che genero nel tempo: ogni colonna rappresenta il campione
%estratto dal Processo di campionamento
Theta=zeros(NVar,T);
 %Punto Iniziale
Theta(1,1)=2; 
Theta(2,1)=2; 
%Campiono
Samples=Fun_MetroPolisHastingsSampler_CW(P,Theta,NVar,T,[1 1]);
zx = sym('zx','real');
zy = sym('zy','real');
Z=[zx;zy];
Zc=[10;10];
F=(Z-Zc).'*(Z-Zc)-16;
ezplot(F,[-20  20 ]); hold on

plot(Samples(1,:),Samples(2,:),'or'); 

