% Addestramento mediante BackPropagation Stocastico
%
clc;clear;close all
format short 
% Definisco la Funzione segno Simbolicamente 
%signum=@(x) x/abs(x)
F=@(x) 1./(1+exp(-x))
%DF=@(x) 1./(1+exp(-x)).^2*exp(-x)
DF=@(x) F(x)*(1-F(x))

%Funzione di Errore
Outpulayer=@(x,a) Neurone([1;Neurone(x,[a(1);a(2);a(3)],'logis');Neurone(x,[a(4);a(5);a(6)],'logis')],[a(7);a(8);a(9)],'logis');
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
a=-1;b=1;
Ho= a + (b-a).*randn(2,3);
H=Ho;
%% Matrice dei Pesi dei Neuroni di Output
% Notazione o(k,j)--> Peso Sinaptico che collega l' input Nascosto xj al
% neurone k di uscita
Oo= a + (b-a).*randn(1,3)
O=Oo;

%Pattern in Ingresso 
%Pattern=[ 1 -1 -1  1;
%          1 -1  1 -1];

Pattern=[ 0 1 1 1;
          0 1 0 0];
      
%Valori di Targhet
t=[0 0 1 1];

%Spazio dei Vettori Augmented
V=Pattern;
V=[ones(1,4);V];
display('Pattern Input')
disp(Pattern)

display('Pesi Hidden')
disp(H)

display('Pesi Output')
disp(O)

H1=H.';
a=[H1(:); O(:)];
for i=1:4
            Err=Errore(V(:,i),a,t(i))
            Out(1,i)=Outpulayer(V(:,i),a);
end
 
display('Uscita')
disp(Out)

display('Errore')
disp(Err)


pause;

E=0;
Tol=1e-2;
J=1;
maxiter=150;
while( (norm(H(1,:))>Tol && norm(H(2,:))>Tol && norm(O)>Tol) && J<=maxiter) 
%for J=1:40 % 5 epoche
    
    display(['Epoca ' num2str(J)])
    teta=.005;
    J=J+1;
    Indici=randperm(4);
    for I=1:4
        %Seleziona il pattern Corrente
        X=V(:,Indici(I));
        
        %Calcolo Netj --> prodotto scalare di Input per i pesi dei Neuroni Hidden
        netj=H*X;
        display('Sinapsi per Neuroni Nascosti')
        disp(netj)
        
        %                                       t
        %Calcolo l' uscita del Layer Hidden: F(W  * X)
        Xj=F(netj);
        display('Uscite Layer Hidden')
        disp(Xj)

        %Aumento lo Spazio delle UScite Hidden(il valore 1 rappresenta il bias per l' unità di Output)
        Xj=[1;Xj];

        %Calcolo Netk --> prodotto scalare di Input per i pesi dei Neuroni Output
        netk=O*Xj;
        display('Sinapsi per Neuroni Output')
        disp(netk)
        
        %Calcolo NetO --> prodotto scalare di Output dei Neuorni Nascosti e dei pesi
        %Calcolo l' uscita del Layer di Output
        Yk=F(netk);
        display('Uscite Layer Output')
        disp(Yk)

        

        %Considero un unita di Uscita qualsiasi


        [mO,nO]=size(O);
        %dE_O_1=zeros(mO,nO);
        for k=1:mO

                %Calcolo la Derivata Parziale di E rispetto al peso di Output o11 : calcolo
                %mediante il METODO CHAIN RULE
                lamda(k)=DF(netk(k))*(Yk(k)-t(k));

                for j=1:nO
                    dE_O(k,j)=lamda(k)*Xj(j);
                end
        end

        %% DERIVATE PARZIALI rispetto ai Neuroni Nascosti
        %Calcolo la Derivata Parziale di E rispetto al peso Hidden H11 : calcolo
        %mediante il METODO CHAIN RULE
        [mH,nH]=size(H);
        dE_h=zeros(mH,nH);

        for j=1:mH
            dE_netj_1=DF(netj(j));
            S=0;
            for k=1:mO
                S=S+lamda(k)*O(k,j);
            end
            lamdaj(j)=S*dE_netj_1;

            for i=1:nH
                dE_h(j,i)=lamdaj(j)*X(i);
            end

        end
       
        %Aggiorna le Matrici dei Pesi
        O=O-teta*dE_O
        H=H-teta*dE_h
        H1=H.';
        
    end
    a=[H1(:); O(:)];
        Out(J,5)=0;
        for i=1:4
            
            Out(J,i)=Outpulayer(V(:,i),a);
            Out(J,5)=Out(J,5)+Errore(V(:,i),a,t(i));
        end
        
end
Out

%Genero un Set di punti in [-1 1] x [-1 1]
[x1 x2]=meshgrid(linspace(-1,1,5));
[m,n]=size(x1);
x=[ones(1,m*n); x1(:).'; x2(:).'];
N=m*n;
for i=1:N
    Z(i)=Outpulayer(x(:,i),a);
    
end


plot(Pattern(1,1:2), Pattern(2,1:2),'or','MarkerSize',8,'MarkerFaceColor','r')
hold on;
plot(Pattern(1,3:4), Pattern(2,3:4),'ob','MarkerSize',8,'MarkerFaceColor','b')

for i=1:N
    DisegnaClusters(x(2,i),x(3,i),Zf(i));
end


