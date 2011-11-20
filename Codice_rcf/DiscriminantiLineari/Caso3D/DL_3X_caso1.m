% DL_3X_caso1.m
%
% Problema di Decisione a 4 classi basato sul criterio MAP con le 
% Verosimiglianze distribuzioni Gaussiana 3-d.
% Le distribuzioni di Verosimiglianza sono distribuzioni gaussiane con la
% stessa Matrice di Varianza - Covarianza -> Sc1=Sc2=Sc3=Sc4 = S = ro^2*I
%
% Procedimento:
%
% 1) Costruisco le Distribuzioni di Verosimiglianza per le 4 classi:
%     P(X | C1), P(X | C2), P(X | C3), P(X | C4) -> Gaussiane 2d aventi la
%     stessa matrice di Varianza - Covarianza
%
% 2) Disegna P(X | C1), P(X | C2), P(X | C3), P(X | C4)
%
% 3) Definisco le Probabilità a priori P(C1), P(C2), P(C3), P(C4)
%
% 4) Costruisco i Discriminanti Lineari
%
% 5) Classifico un insieme di Pattern X
%
% 6) Traccia gli Iperpiani di Separazione - rette in R2


clear;clc;close all;
%Matrice di Varianza Covarianza S=ro*I
ro=5;
S=[ro^2 0   0
   0   ro^2 0 
   0   0    ro^2];
%
% La variabile aleatoria x che utilizzo come caratteristica per la
% classificazione è la (lunghezza,peso).
%
% C1: classe dei Salmoni con Mc1=[50;35]
% 
% C2: classe dei Branzini con Mc2=[70;50]
%
% C3: classe dei Tonni con Mc3=[100;60]
%
% C4: classe dei Merluzzi con Mc4=[65;45;]
%


syms x y z;
X=[x;y;z];
k = menu('Seleziona i vettori delle medie','1','2','3')
 switch k
        case 1
            Mc1=[10;20;15];
            Mc2=[70;50;20];
            Mc3=[80;60;90];
            Mc4=[65;45;50];
        case 2
            Mc1=[10;15;10];
            Mc2=[20;50;20];
            Mc3=[60;30;60];
            Mc4=[50;45;40];
        case 3
            Mc1=[60;15;20];
            Mc2=[20;80;40];
            Mc3=[60;30;50];
            Mc4=[50;45;40];
 end
 

% Costruisco la Distribuzione P(X | C1)
Px_c1=GaussianaMulti(S,Mc1,X);


% Costruisco la Distribuzione P(X | C2)
Px_c2=GaussianaMulti(S,Mc2,X);

% Costruisco la Distribuzione P(X | C3)
Px_c3=GaussianaMulti(S,Mc3,X);

% Costruisco la Distribuzione P(X | C4)
Px_c4=GaussianaMulti(S,Mc4,X);


%legend('P( X | C1)','P( X | C2)','P( X | C3)','P( X | C4)')
 
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
Pc3=0.25;
Pc4=0.25;

%La funzioni discriminanti lineari sono del tipo:
%
% Gi(X)=w(i)'*x+w(i,0)   i=1:4->numero di classi
%
% con:
%
% w(i)=    M(i)/ro^2
%
%
% w(i,0)=  (-1/2*ro^2)*M(i)'*M(i)+log(P(Ci))
%
% 
%
% Il criterio MAP classifica nel seguente modo:
%
% D(x)= max (G1(x),G2(x),G3(x),G4(x))
%          i
%
w=[Mc1 Mc2 Mc3 Mc4].'/ro^2;

wo=(-1/(2*ro^2))*[Mc1.'*Mc1+log(Pc1)
              Mc2.'*Mc2+log(Pc2)
              Mc3.'*Mc3+log(Pc3)
              Mc4.'*Mc4+log(Pc4)];
              
[x,y,z]=meshgrid(linspace(0,150,20));
[m,n,p]=size(x);
x=reshape(x,1,m*n*p);
y=reshape(y,1,m*n*p);
z=reshape(z,1,m*n*p); 

X=[x;y;z];   
size=m*n*p;
G=w*X;
for I=1:size
    G(:,I)=G(:,I)+wo;   
    P1(I)=GaussianaMulti(S,Mc1,X(:,I));
    P2(I)=GaussianaMulti(S,Mc2,X(:,I));
    P3(I)=GaussianaMulti(S,Mc3,X(:,I));
    P4(I)=GaussianaMulti(S,Mc4,X(:,I));
end


x=reshape(x,m,n,p);
y=reshape(y,m,n,p);
z=reshape(z,m,n,p);

P1=reshape(P1,m,n,p);
P2=reshape(P2,m,n,p);
P3=reshape(P3,m,n,p);
P4=reshape(P4,m,n,p);



%Classifica secondo il criterio MAP
[V,I]=max(G);
% Ridistribuisce i dati
DisegnaClusters3D(x(:),y(:),z(:),I);

hold on;

[X,Y,Z]=ellipsoid(Mc1(1),Mc1(2),Mc1(3),15,15,15)
surf(X,Y,Z);

[X,Y,Z]=ellipsoid(Mc2(1),Mc2(2),Mc2(3),15,15,15)
surf(X,Y,Z);

[X,Y,Z]=ellipsoid(Mc3(1),Mc3(2),Mc3(3),15,15,15)
surf(X,Y,Z);

[X,Y,Z]=ellipsoid(Mc4(1),Mc4(2),Mc4(3),15,15,15)
surf(X,Y,Z);

