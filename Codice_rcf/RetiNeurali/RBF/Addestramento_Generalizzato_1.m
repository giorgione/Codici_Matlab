% AddestramentoXor.m
% Problema dello XOR : Addestramento mediante BackPropagation Stocastico 
%                      La Rete è un RBF BIASED
clc;clear;close all
format short 
% Definisco la Funzione Radiale Centrata in t e la sua Derivata 
F=@(x,t,m1,dmax) exp(-(m1/(dmax^2))*(x-t).'*(x-t));


load RBF_data.mat 
plot(P1(1,:),P1(2,:),'or','MarkerSize',5,'MarkerFaceColor','r'); %target 0
hold on
plot(P2(1,:),P2(2,:),'ob','MarkerSize',5,'MarkerFaceColor','b') %target 1

Pattern=[ P1 P2];
 
%Valori di Targhet
Target=[T1 T2];

%Numero Input
Mo=2;
%Numero Neuroni Nascosti
M1=20;
%Numero Neuroni Output
M2=1;

%Numero Training Pattern
N=20;

%Seleziono Casualmente 2 centri dai Dati di Training
I=randperm(20);
t=Pattern(:,I(1:N));
display('Centri Iniziali:')
disp(t)

%Mischio i Pattern
I=randperm(20);
Pattern=Pattern(:,I(1:N));
Target=Target(:,I(1:N));

%t=[P1(:,1) P2(:,2)];
%t=Target;
dmax=norm(t(:,1)-t(:,2));
%calcolo dmax
I=1;
for i=1:M1
    for j=i+1:M1
        D(I)=norm(t(:,1)-t(:,2));
        I=I+1;
    end
end
dmax=max(D)
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


F_RBF=@(x)w.'*[F(x,t(:,1),M1,dmax);F(x,t(:,2),M1,dmax);F(x,t(:,3),M1,dmax);F(x,t(:,4),M1,dmax);...
               F(x,t(:,5),M1,dmax);F(x,t(:,6),M1,dmax);F(x,t(:,7),M1,dmax);F(x,t(:,8),M1,dmax);...
               F(x,t(:,9),M1,dmax);F(x,t(:,10),M1,dmax);F(x,t(:,11),M1,dmax);F(x,t(:,12),M1,dmax);...
               F(x,t(:,13),M1,dmax);F(x,t(:,14),M1,dmax);F(x,t(:,15),M1,dmax);F(x,t(:,16),M1,dmax);...
               F(x,t(:,17),M1,dmax);F(x,t(:,18),M1,dmax);F(x,t(:,19),M1,dmax);F(x,t(:,20),M1,dmax);1] ;
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
[x1 x2]=meshgrid(linspace(-2,2,50));
[m,n]=size(x1);
x=[x1(:).'; x2(:).'];
N=m*n;
for j=1:N
    Z(j)=sign(F_RBF(x(:,j)));
    DisegnaClusters(x(1,j),x(2,j),Z(j))
   
end