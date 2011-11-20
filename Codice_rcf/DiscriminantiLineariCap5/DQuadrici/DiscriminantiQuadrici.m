% DiscriminantiQuadrici.m
%
% Discriminanti Lineari caso Wo = 0
%            t     t
% g(x)=Wo + w x + x  W x  -> FORMA QUADRICA
%
% Per un problema a 2 classi ottengo il seguente comportamento:
%
% 1) g(x) > 0 --> x è in C1 
%
% 2) g(x) < 0 --> x è in C2
%
% 3) g(x) = 0 --> x può appartenere sia a C1 che C2 : non prendo decisione
%                 IPERPIANO DI SEPARAZIONE (quando g(x) è LINEARE)
%
% Procedimento:
%
% 1) Costruisco il classificatore Lineare g(x) a partire dai seguenti
%    valori:
%
%    - w = [2 ;3;1]
%
%   -  W = rand(3)
%
%    - wo = 0 
%
% 2) Disegno il Piano di Separazione determinato da g(x)
%
% 3) Genero un in insieme di Pattern 3-d
%
% 4) Seleziono random 10 Pattern e li classifico
clc;clear;close all
%Vettore dei Pesi
%
w=[2;3;1];
W=[1 2 1;
   2 3 2;
   1 2 1];
wo=0;
%                            2
% Discriminante Lineare_ g: R -> R  
G=@(w,W,x,wo) wo+w.'*x + x.'*W*x;

[x1 x2 x3]=ndgrid(-10:2:10);
[m n p]=size(x1);

D=[x1(:).'; x2(:).'; x3(:).'];

[d,n]=size(D)
for i=1:n
    z(i)=G(w,W,D(:,i),wo);
end

syms x1 x2 x3
x=[x1;x2;x3];

G=wo+w.'*x+x.'*W*x;
Sol=solve(G,'x3');
ezsurf(Sol(1),[-10 10 -10 10]);hold on;
ezsurf(Sol(2),[-10 10 -10 10]);


title('Superfice di Decisione')
xlabel('x1');
ylabel('x2');
zlabel('x3');
hold on;

% Seleziona 10 punti Random e li classifico
indici=round(1+rand(1,20)*n);
for i=1:length(indici)
    I=indici(i);
    if z(I)>0
        %Il punto appartiene alla Regione R1: lato positivo
        plot3(D(1,I),D(2,I),D(3,I),'or','MarkerSize',5,'MarkerFaceColor','r')
    else
        %Il punto appartiene alla Regione R2: lato Negativo
        plot3(D(1,I),D(2,I),D(3,I),'ob','MarkerSize',5,'MarkerFaceColor','b')
    end
end

%Disegno gli assi cartesiani
compas3d(10,0,0,'y')
compas3d(0,10,0,'y')
compas3d(0,0,10,'y')
