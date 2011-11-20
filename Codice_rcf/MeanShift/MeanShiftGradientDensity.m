%Esempio del Mean Shift
%
%Genero dei punti da una Mistura di gaussiane 2D
clc;clear;close all;
mu1 = [2 3];
mu2 = [ 30 10];
SIGMA1 = [10  0.5;
          0.5 3];
SIGMA2 = [10 -0.6;
          -0.6 8 ];
Npoints=100;
Xi =[ mvnrnd(mu1,SIGMA1,Npoints) ;mvnrnd(mu2,SIGMA2,Npoints)];
Npoints=2*Npoints;

%Disegno  il punto iniziale
plot(Xi(:,1),Xi(:,2),'+');hold on;

%Disegno le Mode della distribuzione
plot(mu1(1),mu1(2),'*r');
plot(mu2(1),mu2(2),'*r');

%syms x u s;
%f=exp(-x^2/2);
%df=diff(f,'x')

g=@(x)-0.5*(exp(-(x^2)));
X=[30 20];
plot(X(1),X(2),'or'); hold on
h=0.9;

%Eseguo T iterazioni per vedere dove va il Mean shift
T=50;
t=2;
 
Xp=zeros(T,2);
Xp(1,:)=X;
while t< T
   %Calcolo il Mean Shift Vector
   X=repmat(X,Npoints,1);
   Y=(X-Xi)/h^2;
   Num=0;
   DeNum=0;
   for i=1:Npoints
        y=norm(Y(i,:),2); %*Y(i,:).'
        Num=Num+Xi(i,:)*g(y);
        DeNum=DeNum+g(y);
   end
   MeanShift=Num/DeNum; 
   Xp(t,:)=Xp(t-1,:)-MeanShift;
   X=Xp(t,:);
   t=t+1;
end

plot(Xp(:,1),Xp(:,2),'og-')
