%Campionamento mediante Markov Chains
%
% Transition Function della Catena
%     t    t-1
% P( x  | x   )

%Distribuzione Beta
Beta=@(x,a,b)((x.^(a-1)).*(1-x).^(b-1))./beta(a,b);

%Transition Function utilizzata
TransitionFun=@(xt,xt1) Beta(xt,200*(.9*xt1+0.05),200*(1-0.9*xt1-0.05));

%T: Numero di Iterazioni della catena
t=500;
%K: Numero di Catene
K=3;
Xt=zeros(t,K);
T=1:t;
T=repmat(T',1,K);
Xt(1,:)=unifrnd(0,1,1,K);
syms x;
for t=2:t-1
    
    %Genero dei nuovi campioni da TransitionFun(xt,xt1) dove
    %xt1 è noto
    %xt è incognito --> utilizzo una tecnica di campionamento per estrarre
    %                   campioni dopo che ho sostituito xt1
    
    %Sostituisco Xt1 nella Transition function ed ottengo la funzione di Xt 
    %da campionare 
    a=200*( .9*Xt(t,:)+0.05 );
    b=200*(1-.9*Xt(t,:)-0.05);
    
    Xt(t+1,:)=betarnd(a, b);
    
end

plot(T,Xt,'-')