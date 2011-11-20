%
% Calcolo la Funzione Errore Quadratico MINIMIZZATA dall' Algoritmo di
% BACKWARD PROPAGATION

clc;clear;
% Definisco la Funzione di ATTIVAZIONE
F=@(x) 1./(1+exp(-x))

% Modellizzo la rete costituita da :
%
% 4 Neuroni in Output
% 4 Neuroni Nascosti
% 2 Neuroni in R2
%
% Voglio Calcolare la Funzione Errore

%% INIZIALIZZO LE MATRICI DEI PESI SINAPTICI SIMBOLICHE

% Matrice dei Pesi dei Neuroni Nascosti
% Notazione h(j,i)--> Peso Sinaptico che collega l' input xi al neurone j
syms h11 h12 h21 h22 h31 h32 h41 h42
H=[h11 h12;
   h21 h22;
   h31 h32;
   h41 h42];

%Matrice dei Pesi dei Neuroni di Output
% Notazione o(k,j)--> Peso Sinaptico che collega l' input Nascosto xj al
% neurone k di uscita
syms o11 o12 o13 o14 o21 o22 o23 o24 o31 o32 o33 o34 o41 o42 o43 o44
O=[o11 o12 o13 o14;
   o21 o22 o23 o24;
   o31 o32 o33 o34;
   o41 o42 o43 o44];

%% FASE DI FORWARD: presento un pattern X alla RETE e ne calcolo l' uscita
%% per i vari livelli di OUTPUT
%Pattern in Ingresso
syms x1 x2;
X=[x1;x2];
display('Pattern Input')
disp(X)
E=0;

%Calcolo Netj --> prodotto scalare di Input per i pesi dei Neuroni Hidden
netj=H*X;
display('Sinapsi per Neuroni Nascosti')
disp(netj)
%                                       t
%Calcolo l' uscita del Layer Hidden: F(W  * X)
Xj=F(netj);
display('Uscite Layer Hidden')
disp(Xj)

%Calcolo Netk --> prodotto scalare di Input per i pesi dei Neuroni Output
netk=O*Xj;
display('Sinapsi per Neuroni Output')
disp(netk)
%Calcolo NetO --> prodotto scalare di Output dei Neuorni Nascosti e dei pesi
%Calcolo l' uscita del Layer di Output
Yk=F(netk);
display('Uscite Layer Output')
disp(Yk)

%Costruisco la funzione di Errore sul singolo pattern
t=[1;2;3;4]
E=(1/2)*((Yk-t).'*(Yk-t));

%% Calcolo la Derivata di E rispetto ai PESI SINAPTICI delle unità di
%% USCITA: 
%%
%% dE_O1(k,j): Derivate calcolate con la Tecnica di DErivazione CHAIN RULE
%%
%% dE_O2(k,j): Derivate calcolate mediante jacobian di MATLAB
%%
%% Verifico che dE_O2 = dE_O1

%Considero un unita di Uscita qualsiasi
syms netk_1 
Yk_1=F(netk_1);  %uscita del Neurone di Output

[m,n]=size(O);
dE_O_1=sym(zeros(m,n));
for k=1:m
    
        %Calcolo la Derivata Parziale di E rispetto al peso di Output o11 : calcolo
        %mediante il METODO CHAIN RULE
        dE_netk_1=diff(Yk_1,netk_1)*(Yk(k)-t(k));
        l=dE_netk_1;
        l=subs(l,netk_1,netk(k));
        lamda(k)=l;
        for j=1:n
            dE_O_1(k,j)=lamda(k)*Xj(j);
        end
end

%Calcolo la Derivata Parziale di E rispetto ai pesi di Output : calcolo
%Diretto
dE_O_2=jacobian(E,[o11, o12, o13, o14, o21, o22, o23, o24, o31, o32, o33, o34, o41, o42, o43, o44]);
dE_O_2=reshape(dE_O_2,4,4);
dE_O_2=dE_O_2.';
%Verfico che le Derivate parziali rispetto ai Pesi di Output sono uguali
isequal(dE_O_1,dE_O_2)

%% DERIVATE PARZIALI rispetto ai Neuroni Nascosti con tecnica di
%% DERIVAZIONE DI FUNZIONI COMPOSTE
%Considero un unita Nascosta qualsiasi
syms netj_1 netj_2 netj_3 netj_4

%Considero le 
netJ=[netj_1 ;netj_2; netj_3; netj_4];

% Riconsidero la funzione di Errore
%Calcolo la funzione di Errore come funzione di netj (Attivazioni sinaptiche del layer hidden)
E1=(F(O*F(netJ))-t).^2;
E1=1/2*sum(E1);
E2=subs(E1,{netj_1,netj_2, netj_3, netj_4},{netj(1),netj(2),netj(3),netj(4)});
%Verifico che E ==E2
isequal(E,E2)

[m,n]=size(H);
for i=1:m    
    for j=1:n
        %Calcola la derivata di netj rispetto ai pesi sinaptici che la
        %costituiscono
        dnetj_hij=diff(netj(i),H(i,j));
        
        %Calcolo la derivata di E1 rispetto a netJ
        dE_netj_1=diff(E1,netJ(i));
        dE_netj_1=subs(dE_netj_1,{netj_1,netj_2, netj_3, netj_4},{netj(1),netj(2),netj(3),netj(4)});
        
        %Applico la regola di DERIVAZIONE per FUNZIONI COMPOSTE
        dE_hij_1=dE_netj_1*dnetj_hij;
        
        %Calcolo la Derivata Parziale di E rispetto al peso Hidden H11 : calcolo
        %mediante derivazione diretta
        dE_hij_2=diff(E,H(i,j));
         
        display(['Derivata parziale rispetto ad h(' num2str(i) ',' num2str(j) ')'])
        isequal(factor(dE_hij_1),factor(dE_hij_2))
    end
end

%% DERIVATE PARZIALI rispetto ai Neuroni Nascosti con tecnica di
%% DERIVAZIONE CHAIN RULE

%Calcolo la Derivata Parziale di E rispetto al peso Hidden H11 : calcolo
%mediante il METODO CHAIN RULE
[m,n]=size(H);
dE_h=sym(zeros(m,n));
% 
Xj_1=F(netj_1);  %uscita del Neurone di Output
dE_netj=diff(Xj_1,netj_1);

for j=1:m
    dE_netj_1=subs(dE_netj,netj_1,netj(j));
    sym S;
    S=0;
    for k=1:4
        S=S+lamda(k)*O(k,j);
    end
    lamdaj(j)=S*dE_netj_1;
        
    for i=1:n
        dE_h(j,i)=factor(lamdaj(j)*X(i));
        dE_h_2(j,i)=factor(diff(E,H(j,i)));
        display(['Derivata parziale rispetto ad h' num2str(j) num2str(i)])
        isequal(dE_h(j,i),dE_h_2(j,i))
    end
    
    
end


%J=jacobian(E,[h11, h12, h21, h22, h31, h32, h41, h42,o11, o12, o13, o14, o21, o22, o23, o24, o31, o32, o33, o34, o41, o42, o43, o44,]);
%J=J.';

