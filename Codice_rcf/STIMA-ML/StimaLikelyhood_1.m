% StimaLikelyhood_1.m
%
% Dimostro che la stima a Massima Verosimiglianza della MEDIA è una :
%
% - STIMA NON POLARIZZATA
%
% Procedimento:
%
% 1) Definisco una Distribuzione Gaussiana N(media=5,1) : essa rappresenta la
% distribuzione di cui voglio andare a stimare i parametri (MEDIA)
%
%    FOR i=1:10  
%
%       2) Estraggo Nc=50 campioni dalla Distribuzione.
%
%       3) Stimo la Media attraverso la tecnica a Massima Verosimiglianza:
%          Media(i)
%    END
% 
% 4) Calcolo il Valore Atteso e la Varianza della Variabile Aleatoria Media  
%
% 5) Calcolo il Valore Atteso e la Varianza della Variabile Aleatoria 
%    Errore di Stima della Media   
%
% 6) Disegno la Distribuzione della Variabile Aleatoria stimata e verifico
%    che essa è: 
%                      - una Stima Non Polarizzata
% 
%                      
clear;clc;close all;


sigma=1;
media = 5; 

Media=zeros(1,10);

%Fisso il numero di Campioni per effettuare la Stima
Nc=50;
Nstime=10;
for i=1:Nstime
            
       
        
        %Estraggo i Campioni
        X=media+randn(1,Nc);
        
        %Applico la Stima a MassimaVerosimiglianza di media 
        Media(i)=mean(X);
        
        %Calcolo l' errore- bias
        Errore_M(i)=Media(i)-media;

end

figure;
plot(1:Nstime,abs(Errore_M));
legend('Errore nella Stima della Media')
xlabel('Stima i-esima')
ylabel('Errore di Stima')


%Calcolo il Valore Atteso della Variabile Aleatoria Stimata Media
m=mean(Media);
v=var(Media);

syms x;
F_stimata=GaussianaMulti(v,m,x);
figure;
ezplot(F_stimata)
title('Variabile Aleatoria Stimata mediante stima ML')
hold on;
plot(m,0,'or')
xlabel('Stima della Media')
ylabel('P(Stima della Media)')
legend(['Media= ' num2str(m)])

%Calcolo il Valore Atteso della Variabile Aleatoria Errore di Stima della
%Media e verifico che la Stima della Media è non polarizzata:
%
%      E{Errore_M} = 0   -> Media Nulla
%
%      V{Errore_M} = molto piccolo  -> Varianza molto piccola
m=mean(Errore_M);
v=var(Errore_M);
F_errore=GaussianaMulti(v,m,x);
figure;
ezplot(F_errore)
title('Variabile Aleatoria Errore di Stima');
hold on;
plot(m,0,'or')
xlabel('Errore di Stima')
ylabel('P(Errore di Stima)')
legend(['Media= ' num2str(m)],['Varianza= ' num2str(v)])
