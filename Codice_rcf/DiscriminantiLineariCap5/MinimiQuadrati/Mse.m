%function [J DJ]=Mse(Y,a,b)
% Calcola la funzione dei Minimi Quadrati J(Y,a,b) ed il suo gradiente DJ(Y,a,b)
% dove:
%
% -Y :Vettore dei Pattern (Training Set)
%
% -a :vettore dei pesi (Spazio in cui vado a disegnare la funzione)
%
% -b : vettore dei Margini
%

function [J ,DJ, H]=Mse(Y,a,b)

[Ma Na]=size(a);
if Na>1
    %Devo Disegnare la Superfice dell' Errore
    A=Y*a;
    [m,n]=size(A);
    %restituisce una matrice se a � una matrice
    B=b*ones(1,n);
    M=A-B;
    for i=1:n
         J(i)=norm(M(:,i),2)^2;
    end
    DJ=0;
    H=0;
else
    A=Y*a-b;
    J=norm(A)^2;
    DJ=2*Y.'*(A);
    H=2*Y.'*Y;
end