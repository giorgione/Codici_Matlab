%function J=SquarredError(xo,x)
%
% Funzione Errore Quadratico Medio per un insieme di N osservazioni x:
%
%         N
%         ----                 2
% J(xo)=  \    ||  xo -x(k)  ||
%         /
%         ----
%         k=1
%
% Il Valore di xo che minimizza J è la Media Aritrmetica delle
% Osservazioni x:
%
% argmin J(xo) = mean(x)
%
% parametri in:
%
%   - xo : punto Soluzione
%
%   - x : Insieme delle osservazioni
%
% parametri out:
%
%   - J : funzione Errore Quadratico in xo

%

function J=SquarredError(xo,x)
    [m,n]=size(x);
    A=xo*ones(1,n);
    S=(A-x)*(A-x).';
    J=sum(diag(S)); 