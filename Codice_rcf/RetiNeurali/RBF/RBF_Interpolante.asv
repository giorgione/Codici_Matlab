% function Y=RBF_Interpolante(x,w,X,modello,o)
%
% Valuta la Funzione a Base Radiale Interpolante con centri X  nel punto x.
%
% Parametri:
%
% - x: punto in cui desidero Valutare la Funzione Interpolante
%
% - w: Pesi della Funzione Interpolante.
%
% - X: Centri delle RBF
%
% - modello: Specifica il tipo di RBF che si desidera utilizzare
%
% - o: Varianza del RBF
function Y=RBF_Interpolante(x,w,X,modello,o)

switch modello
    case 1 %Gauss
        T=@(r,o) exp(-(r^2/o^2));
    case 2
        T=@(r,o)((r^2+o^2)^.5)/o^2;
    case 3
        T=@(r,o) o^2/((r^2+o^2)^.5);
    case 4
        T=@(r,o) o^2/(r^2+o^2);
end
n=length(w);
Y=0;
for i=1:n
    Y=Y+w(i)*T(norm(x-X(:,i)),o);
end