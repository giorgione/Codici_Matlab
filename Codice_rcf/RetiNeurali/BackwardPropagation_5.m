% Problema dello XOR : Addestramento mediante BackPropagation Stocastico 
%                      La Rete è NOT BIASED
clc;clear;close all
format short 
% Definisco la Funzione Sigmoidale e la sua Derivata 
F=@(x) 1./(1+exp(-x))
%DF=@(x) 1./(1+exp(-x)).^2*exp(-x)
DF=@(x) F(x)*(1-F(x))

%Funzione di Errore  /norm([a(3);a(4)]
Outpulayer=@(x,a) Neurone( [Neurone(x,[a(1);a(2)]/norm([a(1);a(2)]),'logis');Neurone(x,[a(3);a(4)]/norm([a(3);a(4)]),'logis')],[a(5);a(6)]/norm([a(5);a(6)]),'logis');
Errore=@(x,a,t)1/2*(t-Outpulayer(x,a))^2;


% La rete è costituita da :
%
% 1 Neuroni in Output
% 2 Neuroni Nascosti
% 2 Neuroni in R2
%
% Voglio Calcolare la Funzione Errore

%% Matrice dei Pesi dei Neuroni Nascosti
% Notazione h(j,i)--> Peso Sinaptico che collega l' input xi al neurone j
%
a=-2.38/sqrt(2);
b=2.38/sqrt(2);
%a=-1;
%b=1;
Ho= a + (b-a).*rand(2,2);
Ho=[Ho(1,:)/norm(Ho(1,:));
Ho(2,:)/norm(Ho(2,:))];

H=Ho;

display('Configurazione Iniziale dei Pesi Unità Hidden')
disp(H)

%% Matrice dei Pesi dei Neuroni di Output
% Notazione o(k,j)--> Peso Sinaptico che collega l' input Nascosto xj al
% neurone k di uscita
Oo= (a + (b-a).*rand(1,2));
Oo=Oo/norm(Oo);
O=Oo;
display('Configurazione Iniziale dei Pesi Unità Output:')
disp(O)

pause;
%Pattern in Ingresso 

Pattern=[ -1 1 -1  1;
          -1 1  1 -1];
  
v=var(Pattern,0,2)*[1 1 1 1];
m=mean(Pattern,2)*[1 1 1 1];
Pattern=(Pattern-m)./v;      

%Valori di Targhet
Target=[0 0 1 1];

%Spazio dei Vettori Augmented
V=Pattern;
nP=4;
ao=[reshape(H.',4,1); O(:)];


display('Pattern Input')
disp(Pattern)

display('Pesi Hidden')
disp(H)

display('Pesi Output')
disp(O)

for i=1:nP
            Err=Errore(V(:,i),ao,Target(i))
            Out(1,i)=Outpulayer(V(:,i),ao);
end
 
display('Uscita')
disp(Out)

display('Errore')
disp(Err)


pause;

E=0;
Tol=1e-2;
J=1;
maxiter=100;
alfa=0.6;

[mO,nO]=size(O);
[mH,nH]=size(H);
A=zeros(mH,nH);
B=zeros(mO,nO);
while( (norm(H(1,:))>Tol && norm(H(2,:))>Tol && norm(O)>Tol) && J<=maxiter) 
    
    %display(['Epoca ' num2str(J)])
    teta=0.05;
    J=J+1;
    %Seleziona Randomicamente le coppie (Pattern,Target)
    Indici=randperm(4);
    x=Pattern(:,Indici);
    T=Target(:,Indici);
    
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

            %Aumento lo Spazio delle UScite Hidden(il valore 1 rappresenta il bias per l' unità di Output)
            Xj_Augmented=Xj;

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
            
            %dE_O_1=zeros(mO,nO);
            for k=1:mO

                    %Calcolo la Derivata Parziale di E rispetto al peso di Output o11 : calcolo
                    %mediante il METODO CHAIN RULE
                    lamda(k)=Yk(k)*(1-Yk(k))*(t-Yk(k));

                    for j=1:nO
                        B(k,j)=teta*(1-alfa)*lamda(k)*Xj_Augmented(j)+alfa*B(k,j);
                        dE_O(k,j)=lamda(k)*Xj_Augmented(j);
                    end
            end

            %% DERIVATE PARZIALI rispetto ai Neuroni Nascosti
            %Calcolo la Derivata Parziale di E rispetto al peso Hidden H11 : calcolo
            %mediante il METODO CHAIN RULE
            
            dE_h=zeros(mH,nH);

            for j=1:mH
                %dE_netj_1=DF(netj(j));
                S=0;
                for k=1:mO
                    S=S+lamda(k)*O(k,j);
                end
                lamdaj(j)=Xj_Augmented(j)*(1-Xj_Augmented(j))*S;

                for i=1:nH %numero di Pesi unità nascoste
                    A(j,i)=teta*(1-alfa)*lamdaj(j)*X(i)+alfa*A(j,i);
                    dE_h(j,i)=lamdaj(j)*X(i);
                end
            end
            

            %Aggiorna le Matrici dei Pesi: BackPropagation classico
            %O=O+teta*dE_O;
            %H=H+teta*dE_h;
            %Aggiorna le Matrici dei Pesi: BackPropagation + Momentum
            O=O+B;
            H=H+A;

    end
    
    %Fine Epoca: calcola gli errori        
    a=[reshape(H.',mH*nH,1); O(:)];
    Out(J,5)=0;
    for i=1:4
            
            Out(J,i)=Outpulayer(Pattern(:,i),a);
            Out(J,5)=Out(J,5)+Errore(Pattern(:,i),a,Target(i));
    end
        
end
display('Fine Learning')
pause;
%Genero un Set di punti in [-1 1] x [-1 1]
[x1 x2]=meshgrid(linspace(-1,1,5));
[m,n]=size(x1);
x=[x1(:).'; x2(:).'];
N=m*n;
figure;
for i=1:N
    Z(i)=Outpulayer(x(:,i),a);
    if(Z(i)>0.5)
        DisegnaClusters(x(1,i),x(2,i),1)
    else
        DisegnaClusters(x(1,i),x(2,i),-1)
    end
end


plot(Pattern(1,1:2), Pattern(2,1:2),'or','MarkerSize',8,'MarkerFaceColor','r')
hold on;
plot(Pattern(1,3:4), Pattern(2,3:4),'ob','MarkerSize',8,'MarkerFaceColor','b')

