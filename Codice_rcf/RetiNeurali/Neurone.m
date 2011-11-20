% function [Y]=Neurone(X,W,type)
% Neurone di un Rete Neurale:
%
% Calcola L' uscita di un Pattern passante per il Neurone specificato dai
% Pesi Sinaptici W e dalla funzione di Attivazione type
%
% parametri:
%
% - X: Input Pattern
%
% - W: Pesi Sinaptici
%
% - type: Tipo di funzione di Attivazione 
%
% parametri out:
%
% Y: uscita del Neurone
function [Y]=Neurone(X,W,type)

switch lower(type)
          case {'linear'}
            disp('Method is linear')
            F=@(x)x;
          case 'sign'
            F=@(x) x/abs(x);
          case 'logis'
            F=@(x)1/(1+exp(-x));
          case 'tan'
            F=@(x) tanh(x);
          case 'heaviside'
            F=@(x)heaviside(x);
end
%Calcolo la Dimensione dello Spazio e il numero di Pattern
% [d,n]=size(X);
% for i=1:n
%     %Calcolo  il Valore di Attivazione per il Pattern i-esimo
%     net(i)=W.'*X(:,i);
%     
%     %Calcolo l' Output prodotto dal Neurone per il Pattern i-esimo
%     Y(i)=F(net(i));
% end

%Calcolo  il Valore di Attivazione per il Pattern i-esimo
net=W.'*X;
     
%Calcolo l' Output prodotto dal Neurone per il Pattern i-esimo
Y=F(net);