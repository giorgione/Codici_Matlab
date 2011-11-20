clc;clear;
% Definisco la Funzione segno Simbolicamente 
%signum=@(x) x/abs(x)
%F=@(x) x/abs(x)
%DF=@(x) 1/abs(x)-x/abs(x)^2*signum(x)
F=@(x) 1/(1+exp(-x))
%DF=@(x) 1/(1+exp(-x))^2*exp(-x)
DF=@(x) F(x)*(1-F(x))
%Funzione di Errore
FunErrore=@(x,a) Neurone([1;Neurone(x,[a(1);a(2);a(3)],'sign');Neurone(x,[a(4);a(5);a(6)],'sign')],[a(7);a(8);a(9)],'sign');


% Modellizzo la rete costituita da :
%
% 4 Neuroni in Output
% 4 Neuroni Nascosti
% 2 Neuroni in R2
%
% Voglio Calcolare la Funzione Errore

%% Matrice dei Pesi dei Neuroni Nascosti
% Notazione h(j,i)--> Peso Sinaptico che collega l' input xi al neurone j
%
Ho=randn(2,3);
H=Ho;
%% Matrice dei Pesi dei Neuroni di Output
% Notazione o(k,j)--> Peso Sinaptico che collega l' input Nascosto xj al
% neurone k di uscita
Oo=randn(1,3);
O=Oo;

%Pattern in Ingresso
Pattern=[ 1 -1 -1  1;
          1 -1  1 -1];
      
%Valori di Targhet
t=[-1 -1 1 1];

%Spazio dei Vettori Augmented
V=Pattern;
V=[ones(1,4);V];
display('Pattern Input')
disp(Pattern)

ao=randn(9,1);
a=ao;
for J=1:30 % 5 epoche
    display(['Epoca ' num2str(J)])
    teta=.5/J;

    for I=1:4
        %Seleziona il pattern Corrente
        X=V(:,I);
        tk=t(I);

        Fun=@(a)1/2*(Neurone([1;Neurone(X,[a(1);a(2);a(3)],'sign');Neurone(X,[a(4);a(5);a(6)],'sign')],[a(7);a(8);a(9)],'sign')-tk)^2;
        %Adestramento attraverso il Metodo del Gradiente discendente
        options = optimset('Display','final','TolFun',1e-8,'MaxIter',1,'TolX',1e-8);
        %}{

        
        %Calcolo il minimo
        [a ,Val,output]= fminsearch(Fun,a,options)
        for i=1:4
            FunErrore(V(:,i),a)
        end
    end
end

