%function Pn=Parzen_Window(t,x,hn,Vn,,method)
% Valuta la distribuzione di probabilit� Pn, stima di Parzen dalle
% osservazioni x, sui punti t
%
% parametri in:
%
% - t: intervallo in cui andare a valutare la stima di Parzen
%
% - hn: ampiezza della finestra di Parzen
%
% - Vn: Volume della Regione
%
% -method: specifica il tipo di finestra da utilizzare
%          r-rettangolare   g-gaussiana
%
% parametri out:
%
% - P: la distribuzione di Probabilit� stimata mediante il metodo di Parzen
%
%
function Pn=Parzen_Window(t,x,hn,Vn,method)
n=length(x);
T=length(t);
Pn=zeros(1,T);
 switch lower(method)
          case {'r','R'}
           for j=1:T
                X=t(j);
                for i=1:n
                    Pn(j)=Pn(j)+rettangolo((X-x(i))/hn,1)/Vn;
                end
           end
          case {'g','G'}
              for j=1:T
                X=t(j);
                for i=1:n
                    Pn(j)=Pn(j)+finestra_gaussiana((X-x(i))/hn)/Vn;
                end
                
              end
 end
              Pn=Pn./n;