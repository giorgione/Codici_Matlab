%Verifico le Proprietà della distribuzione Binomiale
clc;close all;clear
a=[0.1 1 2 3];
b=[0.1 1 3 4];

figure;
u=linspace(0,1,100);

for i=1:4
    
        subplot(2,2,i);
        plot(u,BetaDistribution(u, a(i),b(i)));
        xlabel('u')
        ylabel(['Beta(u| ' num2str(a(i)) ',' num2str(b(i)) ]);
     
end

%Considero la Distribuzione BINOMIALE Bin(m|u,N)
%Cerco di stimare la distribuzione di u attraverso stima Bayesiana
%
% Likelyhood= Bin(m|u,N)
% Prior = BetaDistribution(u,a,b)
%
% Post=Prior*Lik
%
% Verifico che Post è ancora una BetaDistribution

%numero di esperimenti
N=10;
%Numero di successi
m=5;
l=N-m;


for i=1:4
    figure(2);
    P= BetaDistribution(u, a(i),b(i)).*DistribuzioneBinomiale(m,N,u);
    subplot(2,2,i);
    plot(u,P);
    xlabel('u')
    ylabel(['Posterior(u| ' num2str(a(i)) ',' num2str(b(i)) ]);
    
    figure(3)
    P1=BetaDistribution(u, m+a(i),l+b(i));
    subplot(2,2,i);
    plot(u,P);
    xlabel('u')
    ylabel(['Posterior(u| ' num2str(a(i)) ',' num2str(b(i)) ]);
end