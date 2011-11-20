% Problema di Decisione a 4 classi basato sul criterio MAP con le 
% Verosimiglianze distribuzioni Gaussiana 3-d.
% Le distribuzioni di Verosimiglianza sono distribuzioni gaussiane con la
%  MATRICE DI VARIANZA-COVARIANZA ARBITRARIA -> Sc1!=Sc2!=Sc3!=Sc4  
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
% 6) Traccia gli Iperpiani di Separazione - curve in R2

clear;clc;close all;
%Matrice di Varianza Covarianza S=ro*I
rx=5;
ry=3;
Sc1=[5^2 0; 
     0   3^2];
 
Sc2=[1^2 0; 
     0   2^2];
 
Sc3=[2^2 0; 
   0   4^2];

Sc4=[4^2 0; 
     0   6^2];
%
% La variabile aleatoria x che utilizzo come caratteristica per la
% classificazione è la (lunghezza,peso,brillantezza).
%
% C1: classe dei Salmoni con Mc1=[50;35] 
% 
% C2: classe dei Branzini con Mc2=[70;50]
%
% C3: classe dei Tonni con Mc3=[100;60]
%
% C4: classe dei Merluzzi con Mc4=[65;45;]
%


syms x y real;
X=[x;y];
k = menu('Seleziona i vettori delle medie','1','2','3');
 switch k
        case 1
            Mc1=[10;20];
            Mc2=[70;50];
            Mc3=[80;60];
            Mc4=[65;45];
            
        case 2
            Mc1=[10;15];
            Mc2=[20;50];
            Mc3=[60;30];
            Mc4=[50;45];
        case 3
            Mc1=[60;15];
            Mc2=[20;80];
            Mc3=[60;30];
            Mc4=[50;45];
 end
M=[Mc1 Mc2 Mc3 Mc4];
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

Pc=[Pc1 Pc2 Pc3 Pc4];


% Costruisco la Distribuzione P(X | C1)
Px_c1=GaussianaMulti(Sc1,Mc1,X);
[F,x,y]=ValutaGaussiana(Sc1,Mc1,-20,100,80);
surf(x,y,F);hold on

% Costruisco la Distribuzione P(X | C2)
Px_c2=GaussianaMulti(Sc2,Mc2,X);
[F,x,y]=ValutaGaussiana(Sc2,Mc2,-20,100,80);
surf(x,y,F)

% Costruisco la Distribuzione P(X | C3)
Px_c3=GaussianaMulti(Sc3,Mc3,X);
[F,x,y]=ValutaGaussiana(Sc3,Mc3,-20,100,80);
surf(x,y,F)

% Costruisco la Distribuzione P(X | C4)
Px_c4=GaussianaMulti(Sc4,Mc4,X);
[F,x,y]=ValutaGaussiana(Sc4,Mc4,-20,100,80);
surf(x,y,F)

legend('P( X | C1)','P( X | C2)','P( X | C3)','P( X | C4)')

%La funzioni discriminanti lineari sono del tipo:
%
% Gi(X)=x'*W(i)*x+w(i)'*x+w(i,0)   i=1:4->numero di classi
%
% con:
%             
% w(i)=    (Sigma(i)^-1)/2
%
% W(i)= -1/2*(Sigma(i)^-1)
%
% w(i,0)=  (-1/2)*M(i)'*(Sigma(i)^-1)*M(i) - 1/2*log(det(Sigma(i))) + log(P(Ci))
%
%
%
% Il criterio MAP classifica nel seguente modo:
%
% D(x)= max (G1(x),G2(x),G3(x),G4(x))
%          i
%
ISc1=inv(Sc1);
ISc2=inv(Sc2);
ISc3=inv(Sc3);
ISc4=inv(Sc4);

w=zeros(2,2,4);
W=zeros(2,2,4);
ISc=zeros(2,2,4);

ISc(:,:,1)=ISc1/-2;
ISc(:,:,2)=ISc2/-2;
ISc(:,:,3)=ISc3/-2;
ISc(:,:,4)=ISc4/-2;

% x=[Mc1(1)+sqrt(Sc1(1,1))*randn(1,20) Mc2(1)+sqrt(Sc2(1,1))*randn(1,20)...
%    Mc3(1)+sqrt(Sc3(1,1))*randn(1,20) Mc3(1)+sqrt(Sc3(1,1))*randn(1,20)];
% y=[Mc1(2)+sqrt(Sc1(2,2))*randn(1,20) Mc2(2)+sqrt(Sc2(2,2))*randn(1,20);...
%    Mc3(2)+sqrt(Sc3(2,2))*randn(1,20) Mc4(2)+sqrt(Sc4(2,2))*randn(1,20)];

[x,y]=meshgrid(linspace(0,150,40), linspace(0,150,40));
[m,n]=size(x);
x=reshape(x,1,m*n);
y=reshape(y,1,m*n);
X=[x;y];  
size=m*n;

for i=1:4

    W(:,:,i)=ISc(:,:,i)/-2;
    
    w(:,i)=ISc(:,:,i)*M(:,i);
    
    wo(i)=(-1/2)*M(:,i).'*ISc(:,:,i)*M(:,i)-(1/2)*log(det(ISc(:,:,i)))+log(Pc(i));
    
    
end



for I=1:size 
    for i=1:4
        G(i,I)=X(:,I).'*W(:,:,i)*X(:,I)*+w(:,i).'*X(:,I)+wo(i);
    end
    P1(I)=GaussianaMulti(Sc1,Mc1,X(:,I));
    P2(I)=GaussianaMulti(Sc2,Mc2,X(:,I));
    P3(I)=GaussianaMulti(Sc3,Mc3,X(:,I));
    P4(I)=GaussianaMulti(Sc4,Mc4,X(:,I));
end


x=reshape(x,m,n);
y=reshape(y,m,n);
P1=reshape(P1,m,n);
P2=reshape(P2,m,n);
P3=reshape(P3,m,n);
P4=reshape(P4,m,n);
figure;

contour(x,y,P1);hold on;plot(Mc1(1),Mc1(2),'ob','MarkerSize',5,'MarkerFaceColor','b');
pause;
contour(x,y,P2);plot(Mc2(1),Mc2(2),'og','MarkerSize',5,'MarkerFaceColor','g');
pause;
contour(x,y,P3);plot(Mc3(1),Mc3(2),'or','MarkerSize',5,'MarkerFaceColor','r');
pause;
contour(x,y,P4);plot(Mc4(1),Mc4(2),'om','MarkerSize',5,'MarkerFaceColor','m');
pause;

%Classifica secondo il criterio MAP
[V,I]=max(G);
% Ridistribuisce i dati
DisegnaClusters(x,y,I);

% P=ginput(1);
% plot(P(1),P(2),'o');
% 
% Punto=[P.' Mc1];
% l=line(Punto(1,:),Punto(2,:));
% set(l,'Color','b','LineWidth',2)
% 
% Punto=[P.' Mc2];
% l=line(Punto(1,:),Punto(2,:));
% set(l,'Color','g','LineWidth',2)
% 
% Punto=[P.' Mc3];
% l=line(Punto(1,:),Punto(2,:));
% set(l,'Color','r','LineWidth',2)
% 
% Punto=[P.' Mc4];
% l=line(Punto(1,:),Punto(2,:));
% set(l,'Color','m','LineWidth',2)

%Calcolo le Superfici di Separazione( curve in questo caso)
% Matlab risolve anche per soluzioni complesse
R12=solve(Px_c1-Px_c2,'x');
l=ezplot(R12(1),'b',[0,150]);set(l,'Color','b','LineWidth',2)
l=ezplot(R12(2),'b',[0,150]);set(l,'Color','b','LineWidth',2)
Punto=[Mc1 Mc2];
l=line(Punto(1,:),Punto(2,:));set(l,'Color','b','LineWidth',2)
axis([0,150,0,150])

pause;

R13=solve(Px_c1-Px_c3,'x');
l=ezplot(R13(1),'r',[0,150]);set(l,'Color','b','LineWidth',2)
l=ezplot(R13(2),'r',[0,150]);set(l,'Color','b','LineWidth',2)
Punto=[Mc1 Mc3];
l=line(Punto(1,:),Punto(2,:));set(l,'Color','b','LineWidth',2)
axis([0,150,0,150])
pause;

R14=solve(Px_c1-Px_c4,'x');
l=ezplot(R14(1),'m',[0,150]);set(l,'Color','b','LineWidth',2)
l=ezplot(R14(2),'m',[0,150]);set(l,'Color','b','LineWidth',2)
Punto=[Mc1 Mc4];
l=line(Punto(1,:),Punto(2,:));set(l,'Color','b','LineWidth',2)
axis([0,150,0,150])
pause;

R23=solve(Px_c2-Px_c3,'x');
l=ezplot(R23(1),'m',[0,150]);set(l,'Color','r','LineWidth',2)
l=ezplot(R23(2),'m',[0,150]);set(l,'Color','b','LineWidth',2)
Punto=[Mc2 Mc3];
l=line(Punto(1,:),Punto(2,:));set(l,'Color','r','LineWidth',2)
axis([0,150,0,150])
pause;

R24=solve(Px_c2-Px_c4,'x');
l=ezplot(R24(1),'m',[0,150]);set(l,'Color','r','LineWidth',2)
l=ezplot(R24(2),'m',[0,150]);set(l,'Color','r','LineWidth',2)
Punto=[Mc2 Mc4];
l=line(Punto(1,:),Punto(2,:));set(l,'Color','r','LineWidth',2)
axis([0,150,0,150])
pause;

R34=solve(Px_c3-Px_c4,'x');
l=ezplot(R34(1),'m',[0,150]);set(l,'Color','m','LineWidth',2)
l=ezplot(R34(2),'m',[0,150]);set(l,'Color','m','LineWidth',2)
Punto=[Mc3 Mc4];
l=line(Punto(1,:),Punto(2,:));set(l,'Color','m','LineWidth',2)
axis([0,150,0,150])

