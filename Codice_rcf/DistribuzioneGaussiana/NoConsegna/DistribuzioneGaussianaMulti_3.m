%Disegna una Gaussiana TREdimensionale , le sue distribuzioni Marginali
% e le distribuzioni condizionate
clc;clear;close all

Ux=2;
Ox=1;

Uy=1;
Oy=.5;

Uz=3;
Oz=1;

M=[Ux;Uy;Uz];
Oxy=0;
Oxz=0;
Oyz=0;
[x,y,z]=ndgrid(-8:.5:8, -8:.5:8, -8:.5:8);

%Matrice di Varianza Covarianza
S=[Ox^2 Oxy Oxz;
   Oxy Oy^2 Oyz;
   Oxz Oyz  Oz^2];
Sinv=inv(S);

[m,n,p]=size(x);
X=[reshape(x,1,m*n*p);reshape(y,1,m*n*p);reshape(z,1,m*n*p)];
F=zeros(m,n,p);

for i=1:m*n*p
    F(i)=GaussianaMulti(S,M,X(:,i));
end
F=reshape(F,m,n,p);

%Disegno P(x) la Distribuzione Marginale della x
Px=sum(sum(F,3),2);
plot(-8:.5:8,Px)
hold on;
plot(Ux,0,'bo','MarkerFaceColor','r');
title('Distribuzione Marginale P(x)');

%Disegno P(y) la Distribuzione Marginale della y
Py=sum(sum(F,3),1);
figure;
plot(-8:.5:8,Py)
hold on;
plot(Uy,0,'bo','MarkerFaceColor','r');
title('Distribuzione Marginale P(y)');

%Disegno P(x) la Distribuzione Marginale della x
Pz=sum(sum(F,1),2);

figure;
plot(-8:.5:8,Pz(:))
hold on;
plot(Uz,0,'bo','MarkerFaceColor','r');
title('Distribuzione Marginale P(z)');

