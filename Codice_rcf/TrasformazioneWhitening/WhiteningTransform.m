%function [M,S,Aw]=WhiteningTransform(M,Sigma)
% Effettua la Trasformata Whitening( reststituisce una nuova Gaussiana con Sigma=I ) 
% per la distribuzione Gaussiana Multivariata avente parametri M ,Sigma:
%
% - Sigma -> Matrice di Varianza-Covarianza
%
% - Media -> Vettore delle Medie
%
% parametri out:
%
% - M -> Vettore delle medie della nuova distribuzione Gaussiana
%
% - S -> Matrice di Varianza-Covarianza della nuova distribuzione Gaussiana
%        Matrice Identita
%
% - Aw -> Matrice che ha generato la trasformazione whitening
function [M,S,Aw]=WhiteningTransform(Media,Sigma)
%calcolo autovalori ed autovettori di Sigma: Lunghezza^2 e Direzioni degli
%assi delle Ellissi
[X,lamda]=eig(Sigma);
%calcolo i vettori Ortonormali di X
[Q,R]=qr(X);
%Genero la Matrice di Trasformazione
Aw=Q*(lamda^(-1/2));

%Trasformo il vettore delle Medie
M=Aw.'*Media;

%Trasformo la Matrice di Varianza-Covarianza -> Matrice Identica
d=length(M);
S=eye(d);