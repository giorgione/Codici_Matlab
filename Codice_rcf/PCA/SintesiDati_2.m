% SintesiDati_2.m
%
% Ragionamento sul significato di MEDIA come Valore rappresentativo per un
% insieme di N osservazioni 3-dimensionali
%
% Procedimento:
%
% 1) Genero un set di N osservazioni 3-D.
%
% 2) Calcolo la Media. 
%
% 3) Calcolo la funzione J - errore quadratico.
%
% 4) Verifico che la Media rappresenta il minimo della Funzione J.
clc;clear ;close all

%Genero un un Set di N dati 3-Dimensionali 
N=10;

a=10;
b=25;
x1=a+(b-a)*rand(1,N);
x2=a+(b-a)*rand(1,N);
x3=a+(b-a)*rand(1,N);

x=[x1;x2;x3];
%Calcolo la media dei Dati: essa rappresenta una Rappresentazione Sintetica
%                           dei dati zero-dimensionale
  
Media=mean(x,2);

%Disegno i dati
plot3(x(1,:),x(2,:),x(3,:),'or','MarkerFaceColor','r');hold on;

%Disegno la Media, intesa come rappresentazione sintetica
plot3(Media(1),Media(2),Media(3),'og','MarkerFaceColor','g');

%La Media rappresenta il valore xo che minimizza la funzione dei minimi
%quadrati:
%         N
%         ----                 2
% J(xo)=  \    ||  xo -x(k)  ||
%         /
%         ----
%         k=1
%
% argmin J(xo) = Media 
%    xo
%ricerco il minimo della funzione J a partire dal punto iniziale x(1):
%ottengo che Sol=Media
[Sol,Val] = fminsearch( @(y) SquarredError(y,x) , [x(1,1);x(2,1);x(2,1)]);

