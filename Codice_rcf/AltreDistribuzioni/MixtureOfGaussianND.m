%Calcolo Distribuzione Mixture of Gaussian
%function L=Mixture(u,Dati,Pw1,Pw2)
%
% parametri:
%
% - X : Matrice delle Osservioni (1 riga-1 Osservazione) m x d
%       m: Numero di Osservazioni
%       d: cardinalità delle osservazioni - num. features

% - M : Matrice delle Medie (n x d) :
%           n: numero di Gaussiane
%           d: cardinalità delle osservazioni - num. features
%
% - S : Matrice 3d delle Covarianze (d x d x n)
%           d: cardinalità delle osservazioni - num. features
%           n: numero di Gaussiane
%
% - Pi : Prior o Mixing coefficient (n x 1)
%
%
% funzione realizzata senza l'uso dei FOR
function P=MixtureOfGaussianND(X,U,S,Pi)
 %Get number of Points --> rows in X
 NumOfGaussian=size(U,1);
 P=zeros(size(X,1),1);
 for i=1:NumOfGaussian
     %Valuta X su tutte le Gaussiane
    P=P+Pi(i)*mvnpdf(X,U(i,:),S(:,:,i));
   
 end 
 