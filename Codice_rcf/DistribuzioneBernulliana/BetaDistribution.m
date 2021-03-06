%%
% Distribuzione Beta: utilizzata come Prior per la stima bayesiana
%                     parametro u nella distribuzione binomiale
%
% 
%function F=Beta(x,a,b)
%Calcola la Distribuzione Beta per la variabile aleatoria u
%
%   P(x|p)= probabilità p che la Variabile Aleatoria X
%             verifichi un successo (x=1) --> p(x=1|p)=p
%             probabilità di sucesso (sempre uguale)
%
%
% PROCESSO DI BERNOULLI:
%
% 1) 2 possibili Risultati mutuamente esclusivi: SUCCESSO - INSUCCESSO
%
% 2) p indica la probabilit� di sucesso per ogni prova: tale valore �
%    costante per ogni prova
%
% 3) Tutti gli Esperimenti - prove sono Indipendenti.
%
% 4) La Variabile Aleatoria X che conta il numero di SUCCESSI su n prove si
%    definisce VARIABILE ALEATORIA BINOMIALE di PARAMETRI p ed n
%
% DISTRIBUZIONE DI PROBABILITA DI BERNULLI:
%
%                   x       1-x
% f(x)= P(X=x) =   p   (1-p)               <-->    P(x|p): parametro della
%                                                          distribuzione
%                 
%
% essa indica:
%             probabilit� che la Variabile Aleatoria X verifichi 1 successo 
%             con p probabilit� di sucesso (sempre uguale)
%
% propriet�:
%
% 1) m = E[ X ] = p    --> Valore MEDIO
% 
%                2
% 2) v = E[ (X-m) ] =  p(1-p)    --> Varianza 


function P=BetaDistribution(u,a,b)
T=@(x)gamma(x);
P=(T(a+b)/(T(a)+T(b)))*(u.^(a-1)).*(1-u).^(b-1);
