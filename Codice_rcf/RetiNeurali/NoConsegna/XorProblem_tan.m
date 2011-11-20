% XorProblem_tan.m
%
% Risolvo il Problema dello Xor in R2 mediante una rete neurale in cui i le
% funzioni di Attivazione per l' HIDDEN LAYER sono di tipo tanh e i cui
% pesi sono gli stessi del Problema risolto mediante percettrone:
%
% IL problema non � risolto poich� la rete non � stata addestrata
%
% Idea:
%
% Ogni Livello della Rete Neurale viene mappato su di un Vettore.
clc;clear;close all
%Dati dello Xor Problem
X=[ 1 -1 -1  1;
    1 -1  1 -1];
n=2;
V=X;
V=[ones(1,2*n);V];

plot(X(1,1:n), X(2,1:n),'or','MarkerSize',8,'MarkerFaceColor','r')
hold on;
plot(X(1,n+1:2*n), X(2,n+1:2*n),'ob','MarkerSize',8,'MarkerFaceColor','b')

HiddenLayer={@(x) Neurone(x,[.5; 1; 1],'tan') ,@(x) Neurone(x,[-1.5; 1; 1],'tan')};

OutputLayer=@(x) Neurone([1;HiddenLayer{1}(x);HiddenLayer{2}(x)],[-1; .7; -.4],'sign');

for i=1:4
    Z(i)=OutputLayer(V(:,i));
end

%Genero un Set di punti in [-1 1] x [-1 1]
[x1 x2]=meshgrid(linspace(-1,1,20));
[m,n]=size(x1);
x=[ones(1,m*n); x1(:).'; x2(:).'];
N=m*n;
for i=1:N
    Z(i)=OutputLayer(x(:,i));
    DisegnaClusters(x(2,i),x(3,i),Z(i))
end

