% Distribuzione di Poisson (o poissoniana) è una distribuzione di probabilità 
% discreta che esprime:
% PROBABILITA  che si verifichino successivamente ed indipendentemente un
% numero di eventi (n) in un dato intervallo di tempo, sapendo che mediamente
% in tale tempo se ne verifica un numero λ: 
%
%                      n         (-λ)
%                     λ       e           
% P(  n | λ ) =  ---------------------   
%                          n!
%
% dove:
%      - n: Numero eventi Indipendenti che si verificano nell'
%           intervallo di Tempo (num. intero)
%
%      -  λ: è il numero medio di eventi che si verifica nell' 
%            intervallo di tempo. (num. intero)
%
% PROPRIETA:
% - E[X]=λ
% - Var[X]=λ
%
% Esempio: 
% Voglio misurare il numero di chiamate ricevute in un call-center in un determinato 
% arco temporale, come una mattinata lavorativa sapendo che mediamente
% ricevo 100 chiamate ogni mattina ---> P(n | 100)
% Questa distribuzione è anche nota come legge degli eventi rari.
%
clc;clear; close all;
PoisonPdf=@(n,lamda) (lamda.^n).*exp(-lamda)./factorial(n)

nt=40;         %punti sui quali valutare la Poison
nlamda=15;     %numero di distr. gamma che voglio generare
t=1:nt;
 
lamda=1:nlamda;
b=linspace(5,20,nlamda);
 
 T=repmat(t,nlamda,1);
 A=repmat(lamda.',1,nt);
 
 P=PoisonPdf(T,A);
 plot(T.',P.','-')
 for i=1:nlamda
     legenda{i}=['lamda=' num2str(lamda(i))];
 end
 legend(legenda)
 
 xlabel( '$n$','Interpreter', 'latex'); 
 ylabel( '$P(n|\lamda)$','Interpreter', 'latex'); 
 

