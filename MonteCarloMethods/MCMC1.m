% Integrazione mediante Montecarlo Methods
% per calcolare il valore Atteso della Distribuzione Beta(x|a,b)
clc;clear;close all
a=3;
b=4;

syms x;
%Calcolo il valore Atteso Analiticamente
Beta=@(x,a,b)((x.^(a-1)).*(1-x).^(b-1))./beta(a,b);

E_x=int(x*Beta(x,a,b),0,1)
disp(['Valore Atteso Analitico: ' num2str(double(E_x))]);
%Campiono la Beta e ne stimo il Valore Atteso Mediante MonteCarlo
%Integration e osservo come al crescere del numero di campioni la
% stima tende al valore esatto <--> AUMENTA L' ACCURATEZZA

N=[100 500 1000 2000];
for i=1:length(N)
    
    X=betarnd(a,b,1,N(i));
    E_xapp(i)=mean(X);
    disp(['Valore Atteso( ' num2str(N(i)) ')='   num2str(E_xapp(i))]);
end


