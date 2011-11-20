clc;clear;
% Definisco la Funzione segno Simbolicamente 
%signum=@(x) x/abs(x)
F=@(x) 1/(1+exp(-x))

% Modellizzo la rete costituita da :
%
% 4 Neuroni in Output
% 4 Neuroni Nascosti
% 2 Neuroni in R2
%
% Voglio Calcolare la Funzione Errore

syms h11 h12 h21 h22 h31 h32 h41 h42
%
% Notazione h(j,i)--> Peso Sinaptico che collega l' input xi al neurone j

%Pesi del Primo Neurone Nascosto
H1=[h11; h12];

%Pesi del secondo Neurone Nascosto
H2=[h21; h22];
%Pesi del terzo Neurone Nascosto
H3=[h31; h32];

%Pesi del 4 Neurone Nascosto
H4=[h41; h42];

%Matrice dei Pesi dei Neuroni Nascosti
H=[h11 h12;
   h21 h22;
   h31 h32;
   h41 h42];

%
% Notazione o(k,j)--> Peso Sinaptico che collega l' input Nascosto xj al
% neurone k di uscita
syms o11 o12 o13 o14 o21 o22 o23 o24 o31 o32 o33 o34 o41 o42 o43 o44

O1=[o11 o12 o13 o14].';
O2=[o21 o22 o23 o24].';
O3=[o31 o32 o33 o34].';
O4=[o41 o42 o43 o44].';

O=[o11 o12 o13 o14;
   o21 o22 o23 o24;
   o31 o32 o33 o34;
   o41 o42 o43 o44];

%Pattern in Ingresso
syms x1 x2;
X=[x1;x2];

%X=[2;1]
E=0;

%Calcolo NetH --> prodotto scalare di Input per i pesi
netH1=H1.'*X;
netH2=H2.'*X;
netH3=H3.'*X;
netH4=H4.'*X;
netH=H*X;
Xj=F(netH);
%                                       t
%Calcolo l' uscita del Layer Hidden: F(W  * X)
Xj1=F(netH1);
Xj2=F(netH2);
Xj3=F(netH3);
Xj4=F(netH4);

%Vettore di Input del Layer Out
Xj=[Xj1;Xj2;Xj3;Xj4];

%Calcolo NetO --> prodotto scalare di Input per i pesi
netO1=O1.'*Xj; 
netO2=O2.'*Xj; 
netO3=O3.'*Xj;
netO4=O4.'*Xj;

%Calcolo l' uscita del Layer di Output
Xk1=F(netO1);
Xk2=F(netO2);
Xk3=F(netO3);
Xk4=F(netO4);

syms net1  net2  net3 net4;
Yk1=F(net1);  %uscita del Neurone 1 di Output
Yk2=F(net2);
Yk3=F(net3);
Yk4=F(net4);


t=rand(1,4);

%Calcolo l' errore sul singolo pattern
E=1/2*[(Xk1-t(1))^2+(Xk2-t(2))^2+(Xk3-t(3))^2+(Xk4-t(4))^2];

E1=1/2*[(Yk1-t(1))^2+(Yk2-t(2))^2+(Yk3-t(3))^2+(Yk4-t(4))^2];

%Calcolo la Derivata Parziale di E rispetto al peso di Output o11 : calcolo
%mediante il METODO CHAIN RULE
dE_net1=diff(Yk1,net1)*(Yk1-t(1));
lamda=dE_net1;
lamda=subs(lamda,net1,netO1);
dE_o11_1=lamda*Xj1;

%Calcolo la Derivata Parziale di E rispetto ai pesi di Output : calcolo
%Diretto
dE_o11_2=diff(E,o11);

isequal(dE_o11_1,dE_o11_2)



J=jacobian(E,[h11, h12, h21, h22, h31, h32, h41, h42,o11, o12, o13, o14, o21, o22, o23, o24, o31, o32, o33, o34, o41, o42, o43, o44,]);
J=J.';

