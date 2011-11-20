%function P=Parzen_Window_Simbolico(X_tr,hn,Vn)
%
% Calcola P la Stima di Parzen Simbolicamente in modo da poterne calcolare l'
% integrale attraverso la funzione int()
%
% parametri in:
%
% - X_tr: campioni estratti dalla distribuzione da stimare
%
% - hn: ampiezza della finestra di Parzen
%
% - Vn: Volume della Regione
%
% parametri out:
%
% - P: la distribuzione di Probabilità stimata mediante il metodo di Parzen
%
function P=Parzen_Window_Simbolico(X_tr,hn,Vn)
syms X real;
P=0;
n=length(X_tr);
for i=1:n
    P=P+finestra_gaussiana((X-X_tr(i))/hn)/Vn;
end
P=P/n;

