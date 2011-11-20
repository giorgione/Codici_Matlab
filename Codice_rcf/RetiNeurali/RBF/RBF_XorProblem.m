% RBF_XorProblem.m
% Risolvo il Problema dello Xor in R2 mediante una Trasormazione non
% Lineare effettuata da 2 Funzioni a BASE RADIALE
%
% Idea:
%
% Ogni Livello della Rete Neurale viene mappato su di un Vettore.
clc;clear;close all
%Dati dello Xor Problem
X=[ 1 0 0 1;
    1 0 1 0];
n=2;


plot(X(1,1:n), X(2,1:n),'or','MarkerSize',8,'MarkerFaceColor','r')
hold on;
plot(X(1,n+1:2*n), X(2,n+1:2*n),'ob','MarkerSize',8,'MarkerFaceColor','b')

t1=[1;1];
t2=[0;0];
F=@(x,t) exp(-((x-t).'*(x-t)));

%Calcolo le trasformazioni non Lineari dei pattern di input
Teta=zeros(4,2);
for i=1:4
    Teta(i,1)=F(X(:,i),t1);
    Teta(i,2)=F(X(:,i),t2);
end


figure;
plot(Teta(:,1),Teta(:,2),'or','MarkerSize',8,'MarkerFaceColor','r');
xlabel('Teta1')
ylabel('Teta2')

syms x1 x2;
X=[x1;x2]
T=[F(X,t1);F(X,t2)];
%Calcolo lo spazio nullo di Teta
figure;
ezsurf(X.'*T)
hold on
plot(Teta(:,1),Teta(:,2),'or','MarkerSize',8,'MarkerFaceColor','r');

