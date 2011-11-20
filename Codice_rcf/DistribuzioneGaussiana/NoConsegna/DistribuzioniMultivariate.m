% Statistica Multivariata
%
% Definisco una funzione di 4 variabili statistiche
% F(x1,x2,x3,x4)
%
% con
%
%       X1 in [1 2]
%
%       X1 in [0 5]
%
%       X1 in [0 12]
%
%       X1 in [2 8]
%

[X1,X2,X3,X4] = ndgrid(1:5, 1:5, 1:5, 1:5);

F=sin(X2.^2).*(X3-2)+exp(X4+X1);
S=sum(sum(sum(sum(F))));

%Faccio in modo che F sia una distribuzione di probabilità
F=F/S;


%Calcolo le distribuzioni Marginali
