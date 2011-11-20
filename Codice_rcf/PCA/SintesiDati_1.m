% SintesiDati_1.m
%
% Ragionamento sul significato di MEDIA come Valore rappresentativo per un
% insieme di osservazioni: esso rappresenta il valore a Minima Distanza 
% (norma-2) da tutte le osservazioni
%
% Procedimento:
%
% 1) Genero un set di N osservazioni 1-D.
%
% 2) Calcolo la Media. 
%
% 3) Calcolo la funzione J - errore quadratico.
%
% 4) Verifico che la Media rappresenta il minimo della Funzione J.
clc;clear ;close all

%Genero un un Set di N dati 1-Dimensionali 
N=10;
a=10;
b=25;
x=a+(b-a)*rand(1,N);

%Calcolo la media dei Dati: essa rappresenta una Rappresentazione Sintetica
%                           dei dati zero-dimensionale

Media=mean(x);

%Disegno i dati
plot(x,ones(1,N),'or','MarkerFaceColor','r');hold on;

%Disegno la Media, intesa come rappresentazione sintetica
plot(Media,1,'og','MarkerFaceColor','g');

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

J=@(xo) sum(norm(xo-x,2).^2)


%ricerco il minimo della funzione J a partire dal punto iniziale x(end):
%ottengo che Sol=Media
[Sol,Val] = fminsearch(J,x(end));

figure;
plot(Sol,Val,'og','MarkerSize',2,'MarkerFaceColor','g');
hold on;
syms t;
F=SquarredError(t,x)
ezplot(F,[a,b]);

