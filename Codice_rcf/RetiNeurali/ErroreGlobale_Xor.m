function F=ErroreGlobale_Xor(X,t,a)
[N d]=size(X);
for i=1:N
    %Seleziono il Pattern i-esimo
    x=X(:,i);
    Outpulayer(i)=Neurone(...
    [1;Neurone(x,[1;a(2);a(3)],'logis');Neurone(x,[1;a(5);a(6)],'logis')]...
    ,[1;a(8);a(9)],'logis');

    Errore(i)=(Outpulayer(i)-t(i))^2;
end

F=sum(Errore);