%function F=DistribuzioneBernulliGen(x,u)
% Calcola la Distribuzione Bernulliana Generalizzata per la variabile
% aleatoria X
%
%   P(x|p)= probabilità p che la Variabile Aleatoria X a n stati
%             verifichi un successo (xi=1) --> p(xi=1|p)=p
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
% 
% 3) Tutti gli Esperimenti - prove sono Indipendenti.
%
% 4) La Variabile Aleatoria X che conta il numero di SUCCESSI su n prove si
%    definisce VARIABILE ALEATORIA BINOMIALE di PARAMETRI p ed n
%
% DISTRIBUZIONE DI PROBABILITA DI BERNULLI GENERALIZZATA:
%
%                        x(k) 
% f(x)= P(X=x) = II  u(k)                  <-->    P(x|p): parametro della
%                                                          distribuzione
%                 
%
% essa indica:
%             probabilit� che la Variabile Aleatoria X verifichi 1 successo 
%             con p probabilit� di sucesso (sempre uguale)
%
% STIMA A MASSIMA LIKELYHOOD:
%
% 1) u(k) = m(k)/N   
% 
function F=DistribuzioneBernulliGen(x,u)

P= u.^x