% DisegnaGaussiana2_sc.m
% 1) Costruisco una Distribuzione Gaussiana Bivariata a variabili
% statisticamente indipendenti
%
% 2) Disegno la Distribuzione Gaussiana Bivariata
%
% 3) Disegno il Grafico delle distribuzioni marginali P(x|y) e P(y|x)
%
% 4) Calcolo gli Auotovalori ed Autovalori della Matrice di Varianza
%    Covarianza e verifico come essi rappresentino rispettivamente la
%    lunghezza e direzione degli assi delle curve di livello della Guassiana
%
clc;clear;close all

Ux=2;
Ox=1;

Uy=1;
Oy=.5;

M=[Ux;Uy];
Oxy=0;
[x y]=meshgrid(-4:.25:8,-4:.25:8);

%Matrice di Varianza Covarianza
S=[Ox^2 Oxy;
   Oxy  Oy^2];

[m,n]=size(x);
X=[reshape(x,1,m*n);reshape(y,1,m*n)];


F=GaussianaMulti(S,M,X);


%Calcola la probabilità condizionata P(y|x=2)
hold on;
[i j]=find(x==2);hold on;
for I=1:length(i)
    z(I)=F(i(I),j(I));
    xI(I)=x(i(I),j(I));
    yI(I)=y(i(I),j(I));
end
plot3(xI,yI,z,'-g','LineWidth',2)

%Calcola la probabilità condizionata P(x|y=1)
[i1 j1]=find(y==1);
for I=1:length(i1)
    z1(I)=F(i1(I),j1(I));
    xI(I)=x(i1(I),j1(I));
    yI(I)=y(i1(I),j1(I));
end
plot3(xI,yI,z1,'-g','LineWidth',2)
xlabel('x');
ylabel('y');

mesh(x,y,F);
Area_Tot=sum(sum(F))



%Calcolo la distribuzione Marginale della X (somma lungo le colonne)
Fx=sum(F,1);
figure;
plot(x(1,:),Fx,'-r','LineWidth',2.2);
hold on;
plot(Ux,0,'bo','MarkerFaceColor','r');
title('Distribuzione Marginale P(x)');
Area=sum(Fx)

%Calcolo la distribuzione Marginale della Y (somma lungo le righe)
Fy=sum(F,2);
figure;
plot(y(:,1),Fy,'-r','LineWidth',2.2);
hold on;
plot(Uy,0,'bo','MarkerFaceColor','r');
title('Distribuzione Marginale P(y)');
Area=sum(Fy)

%Calcolo gli Autovalori ed Autovalori della matrice Sigma
[X,lamda]=eig(S);
figure;
hold on;
a1=sqrt(lamda(1,1));
a2=sqrt(lamda(2,2));

line([X(1,1)+Ux; -X(1,1)+Ux],[X(2,1)+Uy;-X(2,1)+Uy])
line([X(1,2)+Ux; -X(1,2)+Ux],[X(2,2)+Uy;-X(2,2)+Uy])

contour(x,y,F);
xlabel('x');
ylabel('y');



