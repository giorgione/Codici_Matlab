%function outcome=CampionamentoInverso(P,N)
%
% Campionamento Inverso della Distribuzione Discreta P attraverso
% l'estrazione di N campioni
%
% Parametri IN:
% 
%   P: Distribuzione Discreta 1-D 
%   N: Numero campioni che si desidera generare
%
% Parametri OUT:
%
%   outcome: campioni generati
function outcome=CampionamentoInversoDiscreto(P,N,type)
%Caso Discreto
 
    %Genero la Distribuzione Cumulativa di P
    Pcum=cumsum(P);

    %genero N campioni Random a distribuzione uniforme in [0 1]
    Y1=unifrnd(0,1,1,N);
    outcome=zeros(1,N);
    %Campionamento Inverso
    for i=1:N
        j=1;
        while Y1(i) >= Pcum(j)
            j=j+1;          
        end
        outcome(i)=j-1;
    end
 