%Calcolo il Gradiente della Funzione di Errore con BIAS utilizzata dal
%BackwardPropagation SIMBOLICAMENTE E NUMERICAMENTE

clc;clear;
% Definisco la Funzione Sigmoidale Simbolicamente 
F=@(x) 1./(1+exp(-x));
DF=@(x) F(x)*(1-F(x));

% Modellizzo la rete costituita da :
%
% 1 Neuroni in Output
% 2 Neuroni Nascosti
% 2 Neuroni in R2
%
% Voglio Calcolare la Funzione Errore

syms h11 h12 h13 h21 h22 h23
% Matrice dei Pesi dei Neuroni Nascosti
% Notazione h(j,i)--> Peso Sinaptico che collega l' input xi al neurone j
%
H1=[ h11 h12 h13 ;
    h21 h22 h23];
H=rand(2,3)
% Matrice dei Pesi dei Neuroni di Output
% Notazione o(k,j)--> Peso Sinaptico che collega l' input Nascosto xj al
% neurone k di uscita
syms o11 o12 o13 
O1=[ o11 o12 o13];
O=rand(1,3);
%Pattern in Ingresso
syms x1 x2  y1 y2 x3 y3 x4 y4;
X=[1 ;  % augmented
   2;
   2];

%Valori di Target
t=[1];

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
Xj_Aug=[1;Xj];
%Calcolo Netk --> prodotto scalare di Input per i pesi dei Neuroni Output
netk=O*Xj_Aug;
display('Sinapsi per Neuroni Output')
disp(netk)
%Calcolo NetO --> prodotto scalare di Output dei Neuorni Nascosti e dei pesi
%Calcolo l' uscita del Layer di Output
Yk=F(netk);
Yk=Yk.';
display('Uscite Layer Output')
disp(Yk)



%Calcolo Netj --> prodotto scalare di Input per i pesi dei Neuroni Hidden
netj1=H1*X;
display('Sinapsi per Neuroni Nascosti')
disp(netj1)
%                                       t
%Calcolo l' uscita del Layer Hidden: F(W  * X)
Xj1=F(netj1);
display('Uscite Layer Hidden')
disp(Xj1)

Xj1=[1;Xj1];
%Calcolo Netk --> prodotto scalare di Input per i pesi dei Neuroni Output
netk1=O1*Xj1;
display('Sinapsi per Neuroni Output')
disp(netk1)
%Calcolo NetO --> prodotto scalare di Output dei Neuorni Nascosti e dei pesi
%Calcolo l' uscita del Layer di Output
Yk1=F(netk1);
Yk1=Yk1.';
display('Uscite Layer Output')
disp(Yk1)

%Costruisco la funzione di Errore sulla totalità dei PATTERN
E=(1/2)*((Yk1-t).'*(Yk1-t));

%% DERIVATE PARZIALI rispetto ai Neuroni di Output con tecnica di
%% DERIVAZIONE CHAIN RULE

[mO,nO]=size(O)
dE_O_1=sym(zeros(mO,nO));

nP=4;
    %Derivata parziale rispetto al peso o(k,j) per il pattern K-esimo
    for k=1:mO

            %Calcolo la Derivata Parziale di E rispetto al peso di Output o11 : calcolo
            %mediante il METODO CHAIN RULE
            lamda(k)=Yk(k)*(1-Yk(k))*(Yk(k)-t(k));
            for j=1:nO
                dE_O_1(k,j)=Xj_Aug(j)*lamda(k);
            end
    end
    
dE_O_2=jacobian(E,[o11, o12, o13]);
dE_O_2=double(subs(dE_O_2,{o11, o12, o13,h11, h12, h13,h21, h22, h23},{O(1),O(2),O(3)...
    ,H(1,1),H(1,2),H(1,3),H(2,1),H(2,2),H(2,3)}))
dE_O_1=double(dE_O_1)
%Verfico che le Derivate parziali rispetto ai Pesi di Output sono uguali
isequal(dE_O_1,dE_O_2)


%% DERIVATE PARZIALI rispetto ai Neuroni Nascosti con tecnica di
%% DERIVAZIONE CHAIN RULE

%Calcolo la Derivata Parziale di E rispetto al peso Hidden H11 : calcolo
%mediante il METODO CHAIN RULE
[mH,nH]=size(H);
dE_h=sym(zeros(mH,nH));
% 
for j=1:mH+1  %numero Neuroni Nascosti
    sym S;
    S=0;
    for k=1:mO %numero Neuroni Out
        S=S+lamda(k)*O(k,j);
    end
    lamdaj(j)=Xj_Aug(j)*(1-Xj_Aug(j))*S;
        
    for i=1:nH %numero di Pesi unità nascoste
        dE_h(j,i)=double(lamdaj(j)*X(i));
    end
    
    
end
dE_h=double(dE_h(2:end,:))
dE_h_2=jacobian(E,[h11, h12, h13,h21, h22, h23]);
dE_h_2=double(subs(dE_h_2,{o11, o12, o13,h11, h12, h13,h21, h22, h23},{O(1),O(2),O(3)...
    ,H(1,1),H(1,2),H(1,3),H(2,1),H(2,2),H(2,3)}));
dE_h_2=reshape(dE_h_2,nH,mH).'
%Verfico che le Derivate parziali rispetto ai Pesi di Output sono uguali
isequal(dE_h,dE_h_2)

