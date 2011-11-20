% AddestramentoXor.m
% Problema dello XOR : Addestramento mediante BackPropagation Stocastico 
%                      La Rete è un RBF BIASED
clc;clear;close all
format short 
% Definisco la Funzione Radiale Centrata in t e la sua Derivata 
F=@(x,t,m1,dmax) exp(-(m1/(dmax^2))*(x-t).'*(x-t));



% La rete è costituita da :
%
% 1 Neuroni in Output
% 2 Neuroni Nascosti
% 2 Neuroni in R2
%
% Voglio Calcolare la Funzione Errore
%Pattern in Ingresso 

Pattern=[ 1  0  0 1;
          1  1  0 0];
 
%Valori di Targhet
Target=[0 1 0 1];

%Numero Input
Mo=2;
%Numero Neuroni Nascosti
M1=2;
%Numero Neuroni Output
M2=1;

%Numero Training Pattern
N=4;

%Seleziono Casualmente 2 centri dai Dati di Training
I=randperm(N);
t=Pattern(:,I(1:2));
t=[1 0;
   1 0];
dmax=norm(t(:,1)-t(:,2));
%dmax=1;
display('Centri iniziali:')
disp(t)

for i=1:M1
    for j=1:N
       G(j,i)=F(Pattern(:,j),t(:,i),M1,dmax);
    end
end

%Aggiungo La base Radiale costante
G=[G  ones(4,1)];
[U,S,V]=svd(G);
Gpinv1=V*pinv(S)*U.';

%Calcolo la pseudo-inversa di G
Gpinv=pinv(G);
d=Target.';
w=Gpinv*d;

F_RBF=@(x)w.'*[F(x,t(:,1),M1,dmax);F(x,t(:,2),M1,dmax);1] ;
figure;hold on;
for j=1:N
    Z(j)=sign(F_RBF(Pattern(:,j)));
    if Z(j)>0
        plot(Pattern(1,j),Pattern(2,j),'ob','MarkerSize',8,'MarkerFaceColor','b')
    else
        plot(Pattern(1,j),Pattern(2,j),'or','MarkerSize',8,'MarkerFaceColor','r')
    end
end

%Genero un Set di punti in [-1 1] x [-1 1]
[x1 x2]=meshgrid(linspace(0,1,20));
[m,n]=size(x1);
x=[x1(:).'; x2(:).'];
N=m*n;
for j=1:N
    Z(j)=sign(F_RBF(x(:,j)));
    DisegnaClusters(x(1,j),x(2,j),Z(j))
   
end

