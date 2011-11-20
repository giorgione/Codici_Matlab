% Discriminanti Lineari Generalizzati:
% utilizzo gli:
%
%           - augmented weight-vector
%
%           - augmented feature-vector
%       
%     y=g(x)=a(1)+a(2)*x+a(3)*x^2   con a=[1;1;1]
%
% Problema di classificazione 1-D
clc;clear;close all
%Vettore dei Pesi
%
a=[1;1;1];
x=linspace(-10,10,40);
n=length(x);
y=[ones(1,n); x; x.^2]

plot3(y(1,:),y(2,:),y(3,:),'LineWidth',2)
xlabel('x1');
ylabel('x2');
zlabel('x3');

syms x
%Calcolo lo spazio Nullo di y: rappresenta il piano di Separazione
y=[1 x x^2];
N=null(y);
DisegnaPiano([-1;1;0],[1;0;1])

%genero 10 punti random nell' intervallo [-5 5]
x=-5+rand(1,10)*5
G=@(x) 1+x+x.^2;
title('Superfice di Decisione')

hold on;

for I=1:length(x)
    
        if G(x(I))>0
            %Il punto appartiene alla Regione R1: lato positivo
            plot3(1,x(I),x(I)^2,'or','MarkerSize',5,'MarkerFaceColor','r')
        end
        if G(x(I))<0
            %Il punto appartiene alla Regione R2: lato Negativo
            plot3(1,x(I),x(I)^2,'ob','MarkerSize',5,'MarkerFaceColor','b')
        end
        if G(x(I))==0
            %Il punto appartiene alla Regione R2: lato Negativo
            plot3(1,x(I),x(I)^2,'om','MarkerSize',5,'MarkerFaceColor','m')
        end
end


