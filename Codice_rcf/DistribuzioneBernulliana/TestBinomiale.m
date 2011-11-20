%Verifico le Proprietà della distribuzione Binomiale
clc;close all;clear
syms x u
n=10;
%il parametro di probabilità è fisso
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


%% L' Evento complementare di P(m|N,u=0.25) a cosa corrisponde?
% Non avere m successi su N lanci =
%        
%       ------
%       \                                                        M
%  E(w)=   P(n|N,u)= 1 -P(m|N,u)
%       /
%       ------
%       per n!=m
%

[r,c]=size(P1)
for i=1:r
    for j=1:c
        P00(i,j)=sum(P1(i,:))-P1(i,j)
    end
end

figure;
surf(m,N,P00)
xlabel('m')
ylabel('N')
zlabel('P(m!=n|u,N)')

P0=1-P1;
figure;
surf(m,N,P0)
xlabel('m')
ylabel('N')
zlabel('1-P(m|u,N)')

% P0+P1= tutti uno
% P00+P1 = tutti uno
% P00 e P0 sono uguali

%Verifico che P(m|N, u) +P(N-m|N, 1-u) =1 per ogni u
figure;
surf(m,N,P1+P00)
xlabel('m')
ylabel('N')
zlabel('P(x=0|u,N)')
%% EQUIVALENZA   P(N-m|N,1-u) = P(m|N,u)
% P(N-m|N,1-u)  Probabilità di N-m insuccessi su N lanci avendo probabilità di singolo
% successo = u e quindi probabilità di insuccesso =1-u
%
P11=DistribuzioneBinomiale(N-m,N,1-U);
figure;
surf(m,N,P11)
xlabel('m')
ylabel('N')
zlabel('P(x=0|u,N)')

