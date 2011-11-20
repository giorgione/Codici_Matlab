%Calcolo L la Log-Likelyhood function della Mistura Gaussiana
%function L=Mixture(u,Dati,Pw1,Pw2)
%
% parametri:
%
%- o : Vettore delle Varianze
%
% - u : Vettore delle Medie
%
% - Dati : Vettore delle Osservazioni
%
% - Pw1 : Conoscenza a priori
%
% - Pw2 : Conoscenza a priori
%
function L=Mixture(o,u,Dati,Pw1,Pw2)
 
 n=length(Dati);
 
 L=0;
 %Calcolo L la Log-Likelyhood function della Mistura Gaussiana
 for i=1:n
     L=L+Pw1*GaussianaMulti(o(1),u(1),Dati(i)) + Pw2*GaussianaMulti(o(1),u(2),Dati(i)) ;
 end 