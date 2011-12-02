%function F=DistribuzioneBinomiale(x,n,p)
%Calcola la Distribuzione di Binomiale per la variabile aleatoria X
%
%   P(x,n,p)= probabilit� p che la Variabile Aleatoria X
%             verifichi x successi su n prove  --> p(x=1|n,p)=p
%             probabilit� di sucesso (sempre uguale)
%
% Tipologia di Problema: voglio calcolare la probabilit� di ottenere X
% SUCCESSI su n prove (otterr� n-X INSUCCESSI)
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
% DISTRIBUZIONE DI PROBABILITA BINOMIALE o DI BERNULLI:
%
%                 ( n )  x       n-x
% f(x)= P(X=x) =  (   ) p   (1-p)               <-->    P(x|n,p)
%                 ( x )
%
% essa indica:
%             probabilit� che la Variabile Aleatoria X verifichi x successi 
%             su n prove con p probabilit� di sucesso (sempre uguale)
%
% propriet�:
%
% 1) m = E[ X ] = n*p    --> Valore MEDIO
% 
%                2
% 2) v = E[ (X-m) ] =  n*p(1-p)    --> Varianza 



function F=DistribuzioneBinomiale(m,N,u)
%% i valori Neg < 0 non sono ammissibili in quanto identificano M(numero di volte
%  cui ho esito positivo) è  maggiore del Numero totale di prove
Neg=N-m;
[I J]=find(Neg<0);
n=length(I);
%il fattoriale di valori negativi non può essere calcolato per definizione
for i=1:n
    Neg(I(i),J(i))=1;
end

[I J]=find(N<0);
n=length(I);
%il fattoriale di valori negativi non può essere calcolato per definizione
for i=1:n
    N(I(i),J(i))=1;
end

F=factorial(N)./(factorial(m).*factorial(Neg));
F=F.*(u.^m).*(1-u).^(Neg);
%Metto a zero tutte le entrate non coerenti (N<m)
for i=1:n
    F(I(i),J(i))=0;
end



