%ApprossimazioneBernullianaRelativa.m
%Verifico che :
%
% per n-> + inf    e      n*P>10
%
% una distribuzione Binomiale B(k/n,1,P) puo essere approssimata mediante 
% una distribuzione Gaussiana N(u,o^2) avente parametri:
%
% 1) u= P
% 
% 2) o= sqrt(P(1-P)/n) 
%
% Procedimento:
% 
% for I=1:15
%
% 1) Incremento il numero di prove: n=no*(I+1)
%
% 2) Genero La distribuzione di Bernulli P(X,n,0.2)
%
% 3) Calcolo i Valori di Media e Varianza.
%
% 4) Genero la Distribuzione Gaussiana.
%
% 5) Calcolo l' errore di approssimazione in Norma-2
%
% end
% 
clc;close all;clear
P=0.2;
n0=10;
figure;hold on
for I=1:15
    passo=I;
    
    %incremento il numero di Prove
    n=n0*passo
    %numero di successi, valore assumibili:
    %[0-n]--> nessun successo- tutti successi
    x=0:passo:n;
    
    %Calcolo La Distribuzione Binomiale
    ybinomiale=DistribuzioneBernulli(x,n,P)
    y=x/n;
        
    Var=(P*(1-P))/n; %o^2
    Media=P; 
    
    %Calcolo la Distribuzione Gaussiana approssimante
    ygaussiana=exp( (-.5)*( ( (y-Media).^2/Var ) ))./(sqrt(2*pi)*sqrt(Var))
    syms X real
    P=GaussianaMulti(Var,Media,X);
    %Calcolo l' integrale e verifico che esso è 1 (sempre uguale ad 1)
    A=double(int(P,X,-inf,+inf));
    display('Area della Distribuzione di Probabilita N(U=3,O=2)')
    disp(A);
    
    %Calcolo l'errore di approssimazione in norma-2
    errore(I)=norm(ybinomiale-ygaussiana,2);
    
    figure;hold on
    plot(y,ybinomiale,'r');
    plot(y,ygaussiana,'b');
    plot(Media,0,'og','MarkerSize',5,'MarkerFaceColor','g')
    legend(['Distribuzione Binomiale n=' num2str(n) ' P=0.2'],'Distribuzione Gaussiana')
    pause;    
end

figure;
plot(1:15,errore,'r')
title('Errore in norma-2')
xlabel('numero di prove')
ylabel('Errore di Approssimazione')
pause;
close all

% Considerazioni:
%
% Osservando il grafico dell' Errore mi accorgo che n, il numero di prove
% utilizzato è troppo piccolo per ottenere una buona approssimazione della
% distribuzione Bernulliana con una Gaussiana , inoltre non posso aumentare
% n poiche il calcolo del coefficiente binomiale genera overflow
%
% Mi accontento delle proprietà Matematiche della Media e Varianza che mi
% consentono di calcolare tali valori per la Variabile Aleatoria y=x/n.


