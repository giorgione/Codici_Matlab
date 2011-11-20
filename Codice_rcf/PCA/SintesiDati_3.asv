% SintesiDati_3.m
%
% Sintetizzo un insieme di N osservazioni monodimensionali mediante una
% retta passante per la Media di equazione:
%
%               xk= m + a(k)*e  
%
% con:
%
% - m = Media
%
% - e = Direzione della Retta
%          t
% - a(k)= e (x(k)-m) --> proiezione della distanza campione (x(k)-media) sulla retta
%                        (proiezione dello scarto sulla retta)
% Procedimento:
%
% 1) Genero un set di N osservazioni 3-D.
%
% 2) Calcolo la Media. 
%
% 3) Calcolo le proiezioni a(k) dei punti su di una retta aventi direttori 
%    specificati dal vettore e, passante per la Media.
%    Le proiezioni dei punti vengono calcolate minimizzando la distanza ai 
%    minimi quadrati.
%
% 4) Cacolo la miglior direzione e dove proiettare.
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
grid on;

% Defisco la Direzione e della Retta in R3: versore
e=[1;1;1]/norm([1;1;1]);
t=-10:.25:10;
n=length(t)
Retta=[];
for i=1:n
    %Equazione parametrica della Retta
    Retta=[Retta Media+t(i)*e;];
end
%Disegno la Retta, intesa come rappresentazione sintetica
plot3(Retta(1,:),Retta(2,:),Retta(3,:),'LineWidth',2);

%La Media rappresenta il valore xo che minimizza la funzione dei minimi
%quadrati:
%                      xo
%         N            |
%         ----     ------------            2
% J(xo)=  \    ||  Media+a(k)*e -  x(k)  ||
%         /
%         ----
%         k=1
%
% argmin J(xo) = [ a1 a2 .....an] --> i coefficienti ak -->proizioni di
%   xo                                    x(k)-Media sulla Retta
%    
%
%Ricerco il minimo della funzione J a partire da un punto iniziale casuale
%(ottengo che Sol=vettore dei coefficienti):
%
% il parametro incognito da calcolare è il vettore a(k)
% dei coefficienti-proiezioni delle osservazioni e sulla retta passante
% per la Media avente direttore e (fissato)

%
options = optimset('MaxFunEvals',4000,'MaxIter',4000,'TolX',1e-10,'TolFun',1e-8);
[A,Val,exitflag,output] = fminsearch( @(y) SquarredError_Retta(Media,y,e,x) , [rand(1,N)],options);
exitflag
output
%risulta che i coefficienti a(k)
%
% a(k)= e.'*(x(k)-Media)
B=Media*ones(1,N);
a=e.'*(x-B)

%Calcolo le proiezioni delle osservazioni sulla retta specificata dal direttore e
Proj=[];
for i=1:N
    %Proj=[Proj Media+a(i)*e];
    Proj=Media+a(i)*e;

    Punti=[x(:,i) Proj];
    l=line(Punti(1,:),Punti(2,:),Punti(3,:));
    set(l,'Color','r','LineWidth',1)
end
title('Retta - con e fissato')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ricerco la direzione ottima per proiettare i Dati, attraverso la 
% minimizzazione della funzione obbiettivo SquarredError_Retta:
% 
% il parametro incognito da calcolare è il direttore e della retta passante
% per la Media che minimizza la distanza delle osservazioni
[E,Val1,exitflag,output] = fminsearch( @(e) SquarredError_Retta(Media,e.'*(x-Media*ones(1,10))...
,e,x) , [rand(3,1)],options);
exitflag
output


% Ricerco la direzione ottima per proiettare i Dati, attraverso il calcolo
% di Autovalori ed Autovettori della Matrice di Scatter

%Calcolo la Matrice di Scatter
S=MatriceScatter(x,Media);
%Calcolo gli Autovalori ed Autovettori
[X,lamda]=eig(S);

%Calcolo l' autovalore di massimo modulo
[val j]=max(diag(lamda));
e_best=X(:,j)

%Disegno la Retta
e=e_best;
t=-10:.25:20;
n=length(t)
Retta=[];
for i=1:n
    %Equazione parametrica della Retta
    Retta=[Retta Media+t(i)*e;];
end

figure;
%Disegno la Retta, intesa come rappresentazione sintetica
plot3(Retta(1,:),Retta(2,:),Retta(3,:),'b','LineWidth',2);
hold on;

%Disegno i dati
plot3(x(1,:),x(2,:),x(3,:),'or','MarkerFaceColor','r');hold on;

%Disegno la Media, intesa come rappresentazione sintetica
plot3(Media(1),Media(2),Media(3),'og','MarkerFaceColor','g');
grid on;

B=Media*ones(1,N);
a=e.'*(x-B)

%Calcolo le proiezioni delle osservazioni sulla retta specificata dal direttore e
Proj=[];
D=0;
for i=1:N
    %Proj=[Proj Media+a(i)*e];
    Proj=Media+a(i)*e;

    Punti=[x(:,i) Proj];
    l=line(Punti(1,:),Punti(2,:),Punti(3,:));
    set(l,'Color','r','LineWidth',1)
    plot3(Proj(1),Proj(2),Proj(3),'og','MarkerFaceColor','r')
    
    %Calcolo il Valore della Distanza dell' osservazione x(i) dalla
    %Proiezione sulla retta specificata da e: risulta che D=Val1
    D=D+norm(x(:,i)-Proj,2)^2;
end

title('Retta - migliore rappresentazione ai minimi quadrati')