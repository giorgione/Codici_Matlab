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
Xi=Xi.';
Npoints=2*Npoints;

%Disegno  il punto iniziale
plot(Xi(1,:),Xi(2,:),'+');hold on;

%Disegno le Mode della distribuzione
plot(mu1(1),mu1(2),'*r');
plot(mu2(1),mu2(2),'*r');

syms x y;
Z=[x;y]
%f=exp(-x^2/2);
%df=diff(f,'x')

g=@(x)0.5*(exp(-(x^2)));
X=[30;20];
plot(X(1),X(2),'or'); hold on
h=6;

%Eseguo T iterazioni per vedere dove va il Mean shift
T=50;
t=2;
 
Xp=zeros(2,T);
Xp(:,1)=X;
plot(Xp(1,1),Xp(2,1),'og');
while t< T
   %Calcolo il Mean Shift Vector
   Circ= (Z-X).'*(Z-X)-h^2;
   ezplot(Circ,[-10^2 10^2])
   %Cerco i punti che cadono nell'IperSfera Centrata in X
   Sx=[];
   Nx=0; 
   for i=1:Npoints
        if (X-Xi(:,i)).'*(X-Xi(:,i))< h^2
            Sx=[Sx Xi(:,i)];
            Nx=Nx+1;
        end
   end
   
   X=repmat(X,1,Nx);
   %Il mean shift è un Media Locale
   MeanShift=sum(Sx-X,2)/Nx;
   %Fhu=Nx/(N*h^2*pi);
   
   Xp(:,t)=MeanShift+Xp(:,t-1);
   X=Xp(:,t);
   plot(Xp(1,1:t),Xp(2,1:t),'og-');
   t=t+1;
end

%plot(Xp(:,1),Xp(:,2),'og-');
%plot(Xp(1,1),Xp(1,2),'+r')
