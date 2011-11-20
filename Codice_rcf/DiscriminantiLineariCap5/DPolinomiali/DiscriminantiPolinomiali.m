% DiscriminantiPolinomiali.m
%
% Considero il Polinomio di grado 2
%
%  y(x)=a(1)+a(2)*x+a(3)*x^2   con a=[1;1;1]
%
%
%  ed ottengo il discriminante Lineare in y:
% 
%        t   
% g(x)= a  * y(x)  con y(x)=[1 x x^2]
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
%    - a = [1;1;1]
%
% 2) Disegno il Piano di Separazione determinato da g(x)--> PARABOLA
%
% 3) Genero un in insieme di Pattern 1-d
%
% 4) Seleziono random 10 Pattern e li classifico
%
% Problema di classificazione 1-D
clc;clear;close all
%Vettore dei Pesi
a=[-1;1;2];

%genero i punti x
x=linspace(-10,10,40);
n=length(x);

%costruisco la base del polionomio
y=[ones(1,n); x; x.^2]

%Disegno i punti nello Spazio Originale R ed i punti nello Spazio mappato R3
figure;
plot(x,zeros(1,n),'or','MarkerSize',5,'MarkerFaceColor','r');
hold on;
%Mapping da R ad R3
plot3(y(1,:),y(2,:),y(3,:),'LineWidth',2);hold on
plot3(y(1,:),y(2,:),y(3,:),'ob','MarkerSize',5,'MarkerFaceColor','r')
xlabel('x1');
ylabel('x2');
zlabel('x3');
title('Mapping da R ad R^3')
view(3)
figure;

%Mapping da R ad R3
plot3(y(1,:),y(2,:),y(3,:),'LineWidth',2);hold on
xlabel('y1');
ylabel('y2');
zlabel('y3');

%Calcolo lo spazio Nullo di a nello spazio trasformato: rappresenta il
%piano di Separazione
N=null(a.');
DisegnaPiano(N(:,1),N(:,2))
plot3(0,0,0,'og','MarkerSize',7,'MarkerFaceColor','g')
xlabel('y1=1');
ylabel('y2=x');
zlabel('y3=x^2');
%genero 20 punti random nell' intervallo [-5 5]
x=-5+rand(1,20)*10
G=@(x,a) a.'*[1;x;x.^2;]
title('Superfice di Decisione Nello Spazio Trasformato')

hold on;
% Seleziona 20 punti Random e li classifico
%indici=round(1+rand(1,20)*n);
for I=1:length(x)
        G(x(I),a)
        if G(x(I),a)>0
            %Il punto appartiene alla Regione R1: lato positivo
            plot3(1,x(I),x(I)^2,'or','MarkerSize',5,'MarkerFaceColor','r')
        end
        if G(x(I),a)<0
            %Il punto appartiene alla Regione R2: lato Negativo
            plot3(1,x(I),x(I)^2,'ob','MarkerSize',5,'MarkerFaceColor','b')
        end
        if G(x(I),a)==0
            %Il punto appartiene alla Regione R2: lato Negativo
            plot3(1,x(I),x(I)^2,'om','MarkerSize',5,'MarkerFaceColor','m')
        end
end


