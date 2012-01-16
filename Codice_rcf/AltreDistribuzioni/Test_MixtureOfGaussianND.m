%Test  Mixture of Gaussiana Multivariata:
%
%1) Genero la Mistura Multivariata con MixtureOfGaussianND
%
%2) Visualizzo

 %Mistura 2d di 3 Gaussiane
 U=[0 0;10 10; 20 20];
 S=zeros(2,2,3);
 S(:,:,1)=3*diag([1;1]);
 S(:,:,2)=S(:,:,1);
 S(:,:,3)=S(:,:,1);
 Pi=[.3 .4 .3].';
 
 [X Y]=meshgrid(-50:50);
 [m,n]=size(X)
 %Passo le Osservazioni per RIGHE
 x=reshape(X,m*n,1);
 y=reshape(Y,m*n,1);
 XX=[ x y];
 
 ZZ=MixtureOfGaussianND(XX,U,S,Pi);
 Z=reshape(ZZ,m,n);
 surf(X,Y,Z)
 