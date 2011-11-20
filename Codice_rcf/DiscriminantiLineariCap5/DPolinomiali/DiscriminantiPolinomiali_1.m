% DiscriminantiPolinomiali_1.m
%
% Considero il Polinomio di grado 2
%
%  y(x)=a(1)+a(2)*x+a(3)*x^2   con a=[-1;1;2]
%
%
%  ed ottengo il discriminante Lineare in y:
% 
%        t   
% g(x)= a  * y(x)
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
% 1) Costruisco il discriminante polinomiale g(x) a partire dai seguenti
%    valori:
%
%    - a = [-1;1;2]
%
% 2) Disegno il Piano di Separazione determinato da g(x)
%
% 3) Genero un in insieme di Pattern 1-d
%
% 4) Seleziono random 10 Pattern e li classifico
%
% Problema di classificazione 1-D
clc;clear;close all
%Vettore dei Pesi
a=[-1;1;2];

syms x y
%Effettuo mapping non Lineare da R2 a R3
Y=[x;y; x*y];
%Calcolo il discriminante Lineare in R3 (rispetto a Y)
g=a.'*Y;
ezsurf(g)
hold on;
xlabel('y1=x');
ylabel('y2=y');
zlabel('y3=x*y');

%Risolvo g(x)=0 --> Piano di Separazione
H=solve(g);

%Calcolo il piano di Separazione come spazio nullo generato da a
N=null(a.');
DisegnaPiano(N(:,1),N(:,2))

%genero 20 punti random nell' intervallo [-5 5]
x=-5+rand(1,20)*10
G=@(x,a) a.'*[1;x;x.^2;]

title('Superfice di Decisione')

hold on;
% Seleziona 20 punti Random e li classifico
indici=round(1+rand(1,20)*19);
for I=1:length(indici)
        if G(x(I),a)>0
            %Il punto appartiene alla Regione R1: lato positivo
            plot3(1,x(I),x(I)^2,'or','MarkerSize',5,'MarkerFaceColor','r')
        end
        if G(x(I),a)<0
            %Il punto appartiene alla Regione R2: lato Negativo
            plot3(1,x(I),x(I)^2,'ob','MarkerSize',5,'MarkerFaceColor','b')
        end
        if G(x(I),a)==0
            %Il punto appartiene al Piano di Separazione
            plot3(1,x(I),x(I)^2,'om','MarkerSize',5,'MarkerFaceColor','m')
        end
end