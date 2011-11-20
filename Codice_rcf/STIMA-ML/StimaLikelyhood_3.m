%StimaLikelyhood_3.m
%
%Verifico che :
%
% per n-> + inf    
%
% la stima a Massima Verosimiglianza di MEDIA e VARIANZA sono :
%
% - STIMA ASINTOTICAMENTE NON POLARIZZATA
%
% - STIMA CONSISTENTE
%
% Procedimento:
% 
% for I=1:15
%
% 1) Incremento il numero di prove: n=no*(I+1)
%
% 2) Genero n campioni i.i.d provenienti dalla Gaussiana
%
% 3) Stimo i Valori di Media e Varianza.
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

sigma=1;
media = 5; 

syms x;
figure;
F=GaussianaMulti(sigma,media,x);
ezplot(F)
title('Distribuzione Originale');

k=menu('Seleziona la Stima della Varianza','1) Stima Biased','2) Stima Not Biased');
k=k-1;

        
Nstime=300;
X=[];
for i=1:Nstime
    passo=i;
    
    %incremento il numero di Estrazioni
    Nc(i)=n0*passo;
    if(i>1)
        N=Nc(i)-Nc(i-1);
    else
        N=n0;
    end
    %Estraggo i nuovi campioni Campioni
    Xnew=media+randn(1,N);
    
    % L' insieme dei campioni conterrà i campioni vecchi + nuovi campioni
    X=[X Xnew];
    %Applico la Stima a MassimaVerosimiglianza di media 
    Media(i)=mean(X);        
    %Calcolo l' errore- bias
    Errore_M(i)=Media(i)-media;

    %Applico la Stima a MassimaVerosimiglianza di Varianza 
    Varianza(i)=std(X,k);        
    %Calcolo l' errore- bias
    Errore_V(i)=Varianza(i)^2-sigma;

end

figure;
plot(Nc,abs(Errore_M));
legend('Errore nella Stima della Media')
xlabel('Stima i-esima')
ylabel('Errore di Stima')

%Calcolo il Valore Atteso della Variabile Aleatoria Stimata Media
m=mean(Media);
v=std(Media,k);

syms x;
F_stimata=GaussianaMulti(v^2,m,x);
figure;
ezplot(F_stimata)
title('P(Media)-> Variabile Aleatoria Stimata: Media')
hold on;
plot(m,0,'or')
xlabel('Stima della Media')
ylabel('P(Stima della Media)')


%Calcolo il Valore Atteso della Variabile Aleatoria Errore di Stima della
%Media e verifico che la Stima della Media è non polarizzata:
%
%      E{Errore_M} = 0   -> Media Nulla
%
%      V{Errore_M} = molto piccolo  -> Varianza molto piccola
m=mean(Errore_M);
v=std(Errore_M,k);
F_errore=GaussianaMulti(v^2,m,x);
figure;
ezplot(F_errore)
title('P(Errore Media) -> Errore di Stima della Media');
hold on;
plot(m,0,'or')
xlabel('Errore di Stima')
ylabel('P(Errore di Stima)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Ragionamento sulla Varianza
figure;
plot(Nc,abs(Errore_V));
legend('Errore nella Stima della Varianza')
xlabel('Stima i-esima')
ylabel('Errore di Stima ')


%Calcolo il Valore Atteso della Variabile Aleatoria Stimata Varianza
m=mean(Varianza);
v=std(Varianza,k);

syms x;
F_stimata=GaussianaMulti(v^2,m,x);
figure;
ezplot(F_stimata)
title('P(Varianza) -> Variabile Aleatoria Stimata: Varianza ')
hold on;
plot(m,0,'or')
xlabel('Stima della Varianza')
ylabel('P(Stima della Varianza)')

%Calcolo il Valore Atteso della Variabile Aleatoria Errore di Stima della
%Varianza e verifico che la Stima della Media è non polarizzata:
%
%      E{Errore_V} = 0   -> Media Nulla
%
%      V{Errore_V} = molto piccolo  -> Varianza molto piccola
m=mean(Errore_V);
v=std(Errore_V,k);
F_errore=GaussianaMulti(v^2,m,x);
figure;
ezplot(F_errore)
title('P(Errore Varianza) -> Errore di Stima della Varianza');
hold on;
plot(m,0,'or')
xlabel('Errore di Stima della Varianza')
ylabel('P(Errore di Stima della Varianza)')

