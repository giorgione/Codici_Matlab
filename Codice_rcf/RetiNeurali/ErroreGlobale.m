% function F=ErroreGlobale(x,t,a)
%
% Calcolo l' Errore Globale sulla totalità dei Pattern
function F=ErroreGlobale(X,t,a)
[d N]=size(X);
Output=zeros(1,N);
Errore=zeros(1,N);
for I=1:N
    %Seleziono il Pattern i-esimo casualmente
    x=X(:,I);
    Outpulayer(I)=Neurone(...
    [1;Neurone(x,[a(1);a(2);a(3)],'logis');Neurone(x,[a(4);a(5);a(6)],'logis')]...
    ,[a(7);a(8);a(9)],'logis');

    Errore(I)=(Outpulayer(I)-t(I))^2;
end

F=sum(Errore)/2;