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

%
% Notazione o(k,j)--> Peso Sinaptico che collega l' input Nascosto xj al
% neurone k
syms o11 o12 o13 o14 o21 o22 o23 o24 o31 o32 o33 o34 o41 o42 o43 o44

O1=[o11 o12 o13 o14].';
O2=[o21 o22 o23 o24].';
O3=[o31 o32 o33 o34].';
O4=[o41 o42 o43 o44].';

%Pattern in Ingresso
%syms x1 x2;
%X=[x1;x2];

Pattern=rand(2,5);
E=0;
for i=1:5
    
X=Pattern(:,i)
%Calcolo Xi--> prodotto scalare di Input per i pesi
netH1=H1.'*X;
netH2=H2.'*X;
netH3=H3.'*X;
netH4=H4.'*X;

%Calcolo l' uscita del Layer Hidden
Xj1=F(netH1);
Xj2=F(netH2);
Xj3=F(netH3);
Xj4=F(netH4);

%Vettore di Input del Layer Out
Xj=[Xj1;Xj2;Xj3;Xj4]

netO1=O1.'*Xj;
netO2=O2.'*Xj;
netO3=O3.'*Xj;
netO4=O4.'*Xj;

%Calcolo l' uscita del Layer di Output
Xk1=F(netO1);
Xk2=F(netO2);
Xk3=F(netO3);
Xk4=F(netO4);

t=rand(5,4);

%Calcolo l' errore sul singolo pattern
    E=E+1/2*[(Xk1-t(1))^2+(Xk2-t(2))^2+(Xk3-t(3))^2+(Xk4-t(4))^2];
end

% Calcolo il gradiente di di E rispetto il Peso di output o11
%
dE_d011=diff(E,o11)
%pretty(dE_d011)

B=diff(F(netO1))*(Xk1-t(1))*Xj1
%pretty(B)

isequal(dE_d011,B)
