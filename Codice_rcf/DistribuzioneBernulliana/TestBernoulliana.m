%Verifico GRAFICAMENTE le Proprietà della distribuzione Bernoulliana
clc;clear;close all
syms x u

U=linspace(0,1,100);

%P(x=1| u)
P1=DistribuzioneBernulli(1,U);
plot(U,P1)
xlabel('u')
ylabel('P(x=1| u)')

%Verifico che la Bernulliana sia una distribuzione
%Calcolando P(x|u)=P(x=1|u)+P(x=0|u) 

%P(x=0| u)
P0=1-P1;
figure;
plot(U,P0)
xlabel('u')
ylabel('P(x=0|u)')

%Dimostro che P00 e P0 sono equivalenti 
P00=DistribuzioneBernulli(0,U);
figure;
plot(U,P0)
xlabel('u')
ylabel('P(x=0|u)')

%Verifico che P(x|u)= P(x=1| u) +P(x=0| u) =1 per ogni u
figure;
plot(U,P1+P00)
xlabel('u')
ylabel('P(x=1|u)+P(x=0|u)')

%% Calcolo di valore Atteso e Varianza della bernulliana: nn riesce con il
%% calcolo simbolico

%P=DistribuzioneBernulli(x,u);
%Calcolo il Valore Atteso
%Ex=int(x*P,x)
%Calcolo la Varianza
%Varx=int(((x-Ex)^2)*P,x)