%Esempio del Mean Shift CLUSTERING
%
%Genero dei punti da una Mistura di gaussiane 2D
clc;clear;close all;
mu1 = [20 20];
mu2 = [ 30 10];
SIGMA1 = [10  0.5;
          0.5 4];
NPoints=100;

X=mvnrnd(mu1,SIGMA1,NPoints);
Media=mean(X).';
S=cov(X);
x = sym('x','real');
y = sym('y','real');
Z=[x;y];

%Disegno  il punto iniziale
plot(X(:,1),X(:,2),'+');hold on;

%Disegno le Mode della distribuzione
plot(mu1(1),mu1(2),'*r');

[V l]=eig(S);
Ellissi=@(theta)V*3*sqrt(l)*[cos(theta); -sin(theta)]+repmat(Media,1,length(theta));
theta=linspace(0,2*pi,50);
E=Ellissi(theta);
plot(E(1,:),E(2,:));hold on;

plot(Media(1),Media(2),'+b')
