% Problema dello XOR : Addestramento mediante BackPropagation Stocastico
%                      calcolato sull' Errore Globale.
%                      La Rete 2-2-1 � BIASED
clc;clear;close all
format short 
% Definisco la Funzione Sigmoidale Simbolicamente 
F=@(x) 1./(1+exp(-x))
DF=@(x) F(x)*(1-F(x))

%Funzione di Errore
Outpulayer=@(x,a) Neurone([1;Neurone(x,[a(1);a(2);a(3)],'logis');Neurone(x,[a(4);a(5);a(6)],'logis')],[a(7);a(8);a(9)],'logis');
Errore=@(x,a,t)1/2*(t-Outpulayer(x,a))^2;





%% Matrice dei Pesi dei Neuroni Nascosti
% Notazione h(j,i)--> Peso Sinaptico che collega l' input xi al neurone j
%
Ho=[-3.9598 3.5405 1.7735;
   -3.9598 -2.8277 2.772];
 

H=Ho;
display('Configurazione Iniziale dei Pesi Unit� Hidden')
disp(H)

%% Matrice dei Pesi dei Neuroni di Output
% Notazione o(k,j)--> Peso Sinaptico che collega l' input Nascosto xj al
% neurone k di uscita
Oo=[-3.9531 4.1839 3.7222];
O=Oo;
display('Configurazione Iniziale dei Pesi Unit� Output:')
disp(O)

pause;
%Pattern in Ingresso 
Pattern=[ 1 -1 -1  1 2 -2 -2  2 1 -1 -1  1 ;
          1 -1  1 -1 2 -2  2 -2 1 -1  1 -1 ];
      
%Valori di Targhet
Target=[0 0 1 1 0 0 1 1 0 0 1 1];
N=length(Target);

%Spazio dei Vettori Augmented
V=Pattern;
V=[ones(1,N);V];

display('Pattern Input')
disp(Pattern)

display('Pesi Hidden')
disp(H)

display('Pesi Output')
disp(O)

H1=H.';
a=[H1(:); O(:)];
for i=1:N
            Err=Errore(V(:,i),a,Target(i))
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
maxiter=700;
nP=N;

[mO,nO]=size(O);
[mH,nH]=size(H);
A=zeros(mH+1,nH);
B=zeros(mO,nO); 

alfa=0.7;
while( (norm(H(1,:))>Tol && norm(H(2,:))>Tol && norm(O)>Tol) && J<=maxiter) 
%for J=1:40 % 5 epoche
    
    display(['Epoca ' num2str(J)])
    teta=0.025;
    J=J+1;
    Indici=randperm(nP);
    x=V(:,Indici);
    T=Target(:,Indici);
    
    Anew=zeros(mH,nH);
    Bnew=zeros(mO,nO); 
    for I=1:nP
         %Seleziona il pattern I-esimo
         X=x(:,I);
         t=T(I);
        
        %Calcolo Netj --> prodotto scalare di Input per i pesi dei Neuroni Hidden
        netj=H*X;
        %display('Sinapsi per Neuroni Nascosti')
        %disp(netj)
        
        %                                       t
        %Calcolo l' uscita del Layer Hidden: F(W  * X)
        Xj=F(netj);
        %display('Uscite Layer Hidden')
        %disp(Xj)

        %Aumento lo Spazio delle UScite Hidden(il valore 1 rappresenta il bias per l' unit� di Output)
        Xj_Augmented=[1;Xj];

        %Calcolo Netk --> prodotto scalare di Input per i pesi dei Neuroni Output
        netk=O*Xj_Augmented;
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
                lamda(k)=Yk(k)*(1-Yk(k))*(t-Yk(k));

                for j=1:nO
                    dE_O(k,j)=lamda(k)*Xj_Augmented(j);
                    B(k,j)=teta*(1-alfa)*lamda(k)*Xj_Augmented(j)+alfa*B(k,j);
                end
        end

        %% DERIVATE PARZIALI rispetto ai Neuroni Nascosti
        %Calcolo la Derivata Parziale di E rispetto al peso Hidden H11 : calcolo
        %mediante il METODO CHAIN RULE
        
        dE_h=zeros(mH,nH);

        for j=1:mH+1
            %dE_netj_1=DF(netj(j));
            S=0;
            for k=1:mO
                S=S+lamda(k)*O(k,j);
            end
            lamdaj(j)=Xj_Augmented(j)*(1-Xj_Augmented(j))*S;
            
            for i=1:nH %numero di Pesi unit� nascoste
                dE_h(j,i)=lamdaj(j)*X(i);
                A(j,i)=teta*(1-alfa)*lamdaj(j)*X(i)+alfa*A(j,i);
            end
        end
        dE_h=dE_h(2:end,:);
        A1=A(2:end,:);
        
       
        %Aggiorna le Matrici dei Pesi : Algoritmo BackwardPropagation
        %Classico
        %O=O+teta*dE_O;
        %H=H+teta*dE_h;
        %Aggiorna le Matrici dei Pesi: BackPropagation + Momentum
        Bnew=Bnew+B;
        Anew=Anew+A1;
        
    end
    O=O+Bnew;
    H=H+Anew;
    
    %Fine Epoca: calcola gli errori        
    a=[reshape(H.',mH*nH,1); O(:)];
    Out(J,nP+1)=0;
    for i=1:nP
            
            Out(J,i)=Outpulayer(V(:,i),a);
            Out(J,nP+1)=Out(J,nP+1)+Errore(V(:,i),a,Target(i));
    end
        
end

pause;
%Genero un Set di punti in [-1 1] x [-1 1]
[x1 x2]=meshgrid(linspace(-2,2,5));
[m,n]=size(x1);
x=[ones(1,m*n); x1(:).'; x2(:).'];
N=m*n;
figure;
for i=1:N
    Z(i)=Outpulayer(x(:,i),a);
    if(Z(i)>0.5)
        DisegnaClusters(x(2,i),x(3,i),1)
    else
        DisegnaClusters(x(2,i),x(3,i),-1)
    end
end


plot(Pattern(1,1:2), Pattern(2,1:2),'or','MarkerSize',8,'MarkerFaceColor','r')
hold on;
plot(Pattern(1,3:4), Pattern(2,3:4),'ob','MarkerSize',8,'MarkerFaceColor','b')

figure;
plot(2:701,Out(2:end,nP+1));
title('MSE Errore')
xlabel('Epoca')
ylabel('MSE(Epoca)')