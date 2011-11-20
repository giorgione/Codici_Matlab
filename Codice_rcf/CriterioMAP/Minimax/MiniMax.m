% Problema di Decisione a 2 classi basato sul criterio Minimax,
% è nota anche la matrice dei Costi A(2 x 2):
%
% A(i,j): costo associato nell' effettuare l'azione a(i) quando lo stato in
%         natura è C(j).
%
% Suppongo che ci siano solo 2 azioni:
%
%  a1 -> decido che lo stato in natura è C1 
%  a2 -> decido che lo stato in natura è C2
%
% Definisco la matrice 0-1 di perdita, essa associa:
%
% - costo 0 ad una decisione correta    ->     A(1,1)=0=A(2,2)
%
% - costo 1 ad una decisione sbagliata  ->     A(1,2)=1=A(2,1)
%
clear;clc;close all

A=[0 1;
   1 0];
%
%
% C1: classe dei Salmoni
% 
% C2: classe dei Branzini
%
% La variabile aleatoria che utilizzo come caratteristica per la
% classificazione è la lungezza:
%
% Suppongo che le Verosimiglianze siano del tipo:
%
% P(x | C1)= N(50,10)  -> la classe C1 ha una Distribuzione della lunghezza
%                         di tipo Gaussiana con media 50 e Varianza 10

syms x ;
X=[x];
M=0;
%Matrice di Varianza Covarianza
S=[10^2 ];
Px_c1=GaussianaMulti(S,M,X);

%
% P(x | C2)= N(80,5)  -> la classe C1 ha una Distribuzione della lunghezza
%                          di tipo Gaussiana con media 80 e Varianza 5
M=[60];
%Matrice di Varianza Covarianza
S=[5^2 ];
Px_c2=GaussianaMulti(S,M,X);

figure;
h1=ezplot(Px_c1,[0 100]);hold on;
set(h1,'LineWidth',2.2,'Color','r')
h2=ezplot(Px_c2,[0 100]);
set(h2,'LineWidth',2.2,'Color','b')
%
% Supponiamo che la conoscenza a Priori non sia nota (è variabile): occorre
% applicare il criterio minimax per stabilire i valori della probabilità a
% priori.
%
syms a
R= A(1,1)-A(2,2)+(A(2,1)-A(1,1))*int(Px_c1,a,inf)+(A(2,2)-A(1,2))*int(Px_c2,-inf,a)
Pc1=0.7;
Pc2=0.3;

figure;
ezplot(R)
Zero=solve(R)

L=Px_c1/Px_c2;

Rmn=A(1,1)+(A(2,1)-A(1,1))*int(Px_c1,Zero,inf)

syms P1 real

%Calcolo il Separatore delle Regioni ke rende uguali gli Errori di
%Decisione
F=int(Px_c1,a,inf)-int(Px_c2,-inf,a)
Zero=double(solve(F));

%Verifico che F(Zero)=0
double(int(Px_c1,x,Zero,inf))
double(int(Px_c2,x,-inf,Zero))

Rmm=double(A(1,1)+(A(2,1)-A(1,1))*int(Px_c1,Zero,inf))


RP=A(2,2)+(A(1,2)-A(2,2))*int(Px_c2,-inf,Zero)+P1*(...
   (A(1,1)-A(2,2))-(A(2,1)-A(1,1))*int(Px_c1,Zero,inf)...
   -(A(1,2)-A(2,2))*int(Px_c2,-inf,Zero));

   ezplot(RP,[0 100])
   title('P(error)')
   
   RP1=(A(1,1)-A(2,2))-(A(2,1)-A(1,1))*int(Px_c1,Zero,inf)...
   -(A(1,2)-A(2,2))*int(Px_c2,-inf,Zero)

   figure;
   ezplot(RP1,[0 1])
   title('P(error)')
   
   R=A(1,1)*P1*int(Px_c1,-inf,Zero)+A(1,2)*(1-P1)*int(Px_c2,-inf,Zero)+...
       A(2,1)*P1*int(Px_c1,Zero,inf)+A(2,2)*(1-P1)*int(Px_c1,Zero,inf);
   