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

% Voglio Calcolare la Funzione Errore
%Pattern in Ingresso 
display('Inserire i Pattern della classe 1')
P1=[];
figure; 
axis([-2 2 -2 2]);hold on;
for i=1:10
    P=ginput(1);    
    plot(P(1),P(2),'ob','MarkerSize',5,'MarkerFaceColor','b')
    P1=[P1 P.'];
end
T1=zeros(1,10);
pause;

display('Inserire i Pattern della classe 1')
P2=[];
for i=1:10
    P=ginput(1);    
    plot(P(1),P(2),'or','MarkerSize',5,'MarkerFaceColor','r')
    P2=[P2 P.'];
end
T2=ones(1,10);


Pattern=[ P1 P2];
 
%Valori di Targhet
Target=[T1 T2];

%Numero Input
Mo=2;
%Numero Neuroni Nascosti
M1=2;
%Numero Neuroni Output
M2=1;

%Numero Training Pattern
N=20;

%Seleziono Casualmente 2 centri dai Dati di Training
I=randperm(N);
%t=Pattern(:,I(1:2));
t=[P1(:,1) P2(:,2)];

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
G=[G  ones(N,1)];
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
[x1 x2]=meshgrid(linspace(-2,2,20));
[m,n]=size(x1);
x=[x1(:).'; x2(:).'];
N=m*n;
for j=1:N
    Z(j)=sign(F_RBF(x(:,j)));
    DisegnaClusters(x(1,j),x(2,j),Z(j))
   
end