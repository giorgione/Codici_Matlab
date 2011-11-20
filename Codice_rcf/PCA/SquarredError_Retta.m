%function J=Squarred_Retta(xo,x)
%
% Funzione Errore Quadratico Medio per un insieme di N osservazioni x:
%
%         N
%         ----                             2
% J(xo)=  \    ||  Media+a(k)*e -  x(k)  ||
%         /
%         ----
%         k=1
%
% Il Valore di xo che minimizza J è la Media Aritrmetica delle
% Osservazioni x:
%
% argmin J(xo) = mean(x)

function J=SquarredError_Retta(Media,a,e,x)
    [m,n]=size(x);
    A=zeros(m,n);
    for i=1:n
        A(:,i)=Media+a(i)*(e/norm(e,2));
    end
    S=(A-x)*(A-x).';
    J=sum(diag(S)); 