% Problema di Decisione a 4 classi basato sul criterio MAP con le 
% Verosimiglianze distribuzioni Gaussiana 3-d.
% Le distribuzioni di Verosimiglianza sono distribuzioni gaussiane con la
% stessa Matrice di Varianza - Covarianza -> Sc1=Sc2=Sc3=Sc4 = S = ro^2*I
%
clear;clc;
%Matrice di Varianza Covarianza S=ro*I
ro=5;
S=[ro^2 0   0;
   0   ro^2  0
   0    0   ro^2];
%
% La variabile aleatoria x che utilizzo come caratteristica per la
% classificazione è la (lunghezza,peso,brillantezza).
%
% C1: classe dei Salmoni con Mc1=[50;35;100]
% 
% C2: classe dei Branzini con Mc2=[70;50;120]
%
% C3: classe dei Tonni con Mc3=[100;60;80]
%
% C4: classe dei Merluzzi con Mc4=[65;45;110]
%

 

syms x y z;
X=[x;y;z];
Mc1=[10;20;20];
Mc2=[70;50;120];
Mc3=[80;60;80];
Mc4=[65;45;110];

Px_c1=GaussianaMulti(S,Mc1,X);

Px_c2=GaussianaMulti(S,Mc2,X);

Px_c3=GaussianaMulti(S,Mc3,X);

Px_c4=GaussianaMulti(S,Mc4,X);


%
% Supponiamo che la conoscenza a Priori sia:
%
% " Nel Mar Baltico il 20% della popolazione dei pesci è Salmone ed il 30%
% è Branzino, il 15% è Tonno, il 35% è Merluzzo "
%
% che si traduce in:
%
% P(C1)=2/10    P(C2)=3/10  P(C3)=15/100     P(C2)=35/100
Pc1=0.25;
Pc2=0.25;
Pc3=0.15;
Pc4=0.35;

%La funzioni discriminanti lineari sono del tipo:
%
% Gi(X)=w(i)'*x+w(i,0)   i=1:4->numero di classi
%
% con:
%
% w(i)=    M(i)/ro^2
%
%
% w(i,0)=  (-1/ro^2)*M(i)'*M(i)+log(P(Ci))
%
%
%
% Il criterio MAP classifica nel seguente modo:
%
% D(x)= max (G1(x),G2(x),G3(x),G4(x))
%          i
%
w=[Mc1 Mc2 Mc3 Mc4].'/ro^2;

wo=(-1/ro^2)*[Mc1.'*Mc1+log(Pc1)
              Mc2.'*Mc2+log(Pc2)
              Mc3.'*Mc3+log(Pc3)
              Mc4.'*Mc4+log(Pc4)];
              
[x,y,z]=ndgrid(linspace(0,100,20), linspace(0,100,20), linspace(0,100,20));
[m,n,p]=size(x);
x=reshape(x,1,m*n*p);
y=reshape(y,1,m*n*p);
z=reshape(z,1,m*n*p)
X=[x;y;z];  
size=m*n*p;
G=w*X;
for I=1:size
    G(:,I)=G(:,I)+wo;
end
    


%Classifica secondo il criterio MAP
[V,I]=max(G)
% Ridistribuisce i dati
%DisegnaClusters(x,y,z,I);

%Calcolo l'iperpiano
H=G(1,:)-G(3,:)
