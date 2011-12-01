%Verifico le Proprietà della distribuzione Binomiale
% DISTRIBUZIONE DI PROBABILITA BINOMIALE:
%
%                 ( n )  x       n-x
% f(x)= P(X=x) =  (   ) p   (1-p)               <-->    P(x|n,p)
%                 ( x )
%
% essa indica:
%             probabilit� che la Variabile Aleatoria BERNOULLIANA X verifichi x successi 
%             su n prove con p probabilit� di sucesso (sempre uguale)
%
% propriet�:
%
% 1) m = E[ X ] = n*p    --> Valore MEDIO
% 
%                2
% 2) v = E[ (X-m) ] =  n*p(1-p)    --> Varianza 
clc;close all;clear
syms x u
n=20;
%il parametro di probabilità è fisso : PROBABILITA DI SUCCESSO
U=0.25;
N=1:n;
m=0:n;
[m,N]=meshgrid(m,N);
%N=10*ones(n+1,n+1);

% Probabilità di m successi su N lanci avendo probabilità di singolo
% successo = u
%P(m|N=10,u=0.25)
P1=DistribuzioneBinomiale(m,N,U);
surf(m,N,P1)
xlabel('m')
ylabel('N')
zlabel('P(x=1 per m volte|u,N)')

%Verifico che P1 sia una distribuzione
%sommo le righe di P1: ogni riga corrisponde ad un differente valore di m
sum(P1,2)


%% L' Evento complementare Pc di P(x=m|N,u=0.25) a cosa corrisponde?
% Non avere m successi su N lanci P(x!=m|N,u=0.25)
%        

[r,c]=size(P1)
for i=1:r
    for j=1:c
        Pc(i,j)=1-P1(i,j);
    end
end

figure;
surf(m,N,Pc)
xlabel('m: numero Insuccessi')
ylabel('N: numero prove effettuate')
zlabel('P(m!=n|u,N)')

P0=1-P1;
figure;
surf(m,N,P0)
xlabel('m: numero Insuccessi')
ylabel('N: numero prove effettuate')
zlabel('1-P(m|u,N)')

% P0+P1= tutti uno
% Pc+P1 = tutti uno
% Pc e P0 sono uguali

%Verifico che P(m|N, u) +P(n!=m|N, 1-u) =1 per ogni u
figure;
surf(m,N,P1+Pc)
xlabel('m: numero successi')
ylabel('N: numero prove effettuate')
zlabel('P(x=0|u,N)')
%% EQUIVALENZA   P(N-m|N,1-u) = P(m|N,u)
% P(N-m|N,1-u)  Probabilità di N-m insuccessi su N lanci avendo probabilità di singolo
% successo = u e quindi probabilità di insuccesso =1-u
%
P11=DistribuzioneBinomiale(N-m,N,1-U);
figure;
surf(m,N,P11)
xlabel('m: numero successi')
ylabel('N: numero prove effettuate')
zlabel('P(x=0|u,N)')

