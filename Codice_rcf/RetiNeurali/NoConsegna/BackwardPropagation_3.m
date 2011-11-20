% Addestramento mediante BackPropagation Batch
%
clc;clear;
% Definisco la Funzione segno Simbolicamente 
%signum=@(x) x/abs(x)
F=@(x) 1/(1+exp(-x))
DF=@(x) 1/(1+exp(-x))^2*exp(-x)

%Funzione di Errore
Outpulayer=@(x,a) Neurone([1;Neurone(x,[a(1);a(2);a(3)],'sign');Neurone(x,[a(4);a(5);a(6)],'sign')],[a(7);a(8);a(9)],'sign');
Errore=@(x,a,t)1/2*(Outpulayer(x,a)-t)^2;

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
Ho=rand(2,3);
H=Ho;
%% Matrice dei Pesi dei Neuroni di Output
% Notazione o(k,j)--> Peso Sinaptico che collega l' input Nascosto xj al
% neurone k di uscita
Oo=rand(1,3);
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


E=0;
Tol=1e-5;
J=1;
maxiter=80;
while( (norm(H)>Tol || norm(O)>Tol) && J<=maxiter) 
%for J=1:40 % 5 epoche
    
    display(['Epoca ' num2str(J)])
    teta=1/J;
    J=J+1;
    Indici=randperm(4);
    
    [mO,nO]=size(O);
    dE_O=zeros(mO,nO);
   
    [mH,nH]=size(H);
    dE_h=zeros(mH,nH);
    
    for I=1:4
        %Seleziona il pattern Corrente
        X=V(:,Indici(I));
        
        %Calcolo Netj --> prodotto scalare di Input per i pesi dei Neuroni Hidden
        netj=H*X;
        %display('Sinapsi per Neuroni Nascosti')
        %disp(netj)
        
        %                                       t
        %Calcolo l' uscita del Layer Hidden: F(W  * X)
        Xj=F(netj);
        %display('Uscite Layer Hidden')
        %disp(Xj)

        %Aumento lo Spazio delle UScite Hidden(il valore 1 rappresenta il bias per l' unità di Output)
        Xj=[1;Xj];

        %Calcolo Netk --> prodotto scalare di Input per i pesi dei Neuroni Output
        netk=O*Xj;
        %display('Sinapsi per Neuroni Output')
        %disp(netk)
        
        %Calcolo NetO --> prodotto scalare di Output dei Neuorni Nascosti e dei pesi
        %Calcolo l' uscita del Layer di Output
        Yk=F(netk);
        %display('Uscite Layer Output')
        %disp(Yk)

        

        %Considero un unita di Uscita qualsiasi
        for k=1:mO

                %Calcolo la Derivata Parziale di E rispetto al peso di Output o11 : calcolo
                %mediante il METODO CHAIN RULE
                lamda(k)=DF(netk(k))*(Yk(k)-t(k));

                for j=1:nO
                    dE_O(k,j)=dE_O(k,j)+lamda(k)*Xj(j);
                end
        end

        

        %% DERIVATE PARZIALI rispetto ai Neuroni Nascosti
        %Calcolo la Derivata Parziale di E rispetto al peso Hidden H11 : calcolo
        %mediante il METODO CHAIN RULE
        

        for j=1:mH
            dE_netj_1=DF(netj(j));
            S=0;
            for k=1:mO
                S=S+lamda(k)*O(k,j);
            end
            lamdaj(j)=S*dE_netj_1;

            for i=1:nH
                dE_h(j,i)=dE_h(j,i)+lamdaj(j)*X(i);
            end

        end
        
    end
    
    %fine Epoca
    %Aggiorna le Matrici dei Pesi
    O=O-teta*dE_O;
    H=H+teta*dE_h;
    H1=H.';

    display('Classificazione')
    for i=1:4
        %Errore(V(:,i),[H1(:); O(:)],t(i))
        Outpulayer(V(:,i),[H1(:); O(:)])
    end
    
end
