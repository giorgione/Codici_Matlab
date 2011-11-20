% Discriminanti Lineari Generalizzati:
% utilizzo gli:
%
%           - augmented weight-vector
%
%           - augmented feature-vector
%       t
% g(x)=W x + Wo
clc;clear;close all

%Vettore dei Pesi
%
wo=8;
w=[wo;2;3;1];
display('Vettore w Augmented')
disp(w)


%                            2
% Discriminante Lineare_ g: R -> R  
G=@(w,x) w.'*x

[x1 x2 x3]=ndgrid(-10:10);
[m n p]=size(x1);

D=[ones(1,m*n*p); x1(:).'; x2(:).'; x3(:).'];
z=G(w,D);

syms x1 x2 x3
x=[1;x1;x2;x3];

G=w.'*x
Sol=solve(G,'x3');
ezsurf(Sol);

title(['Superfice di Decisione : wo=' num2str(wo)])
xlabel('x1');
ylabel('x2');
zlabel('x3');
hold on;
% Seleziona 10 punti Random e li classifico
indici=round(1+rand(1,10)*m*n*p);
for i=1:length(indici)
    I=indici(i);
    
        if z(I)>0
            %Il punto appartiene alla Regione R1: lato positivo
            plot3(D(1,I),D(2,I),D(3,I),'or','MarkerSize',5,'MarkerFaceColor','r')
        end
        if z(I)<0
            %Il punto appartiene alla Regione R2: lato Negativo
            plot3(D(1,I),D(2,I),D(3,I),'ob','MarkerSize',5,'MarkerFaceColor','b')
        end
        if z(I)==0
            %Il punto appartiene alla Regione R2: lato Negativo
            plot3(D(1,I),D(2,I),D(3,I),'om','MarkerSize',5,'MarkerFaceColor','m')
        end
end


Compas3d(10,0,0,'y')
Compas3d(0,10,0,'y')
Compas3d(0,0,10,'y')
Compas3d(w(1),w(2),w(3),'m')

%Disegno 2 punti X1 e X2 appartenenti all'Iperpiano H di separazione
A=subs(Sol,{x1,x2},{1,1});
X1=[1;1;A];
plot3(X1(1),X1(2),X1(3),'om','MarkerSize',5,'MarkerFaceColor','m')
B=subs(Sol,{x1,x2},{2,3});
X2=[2;3;B];
plot3(X2(1),X2(2),X2(3),'om','MarkerSize',5,'MarkerFaceColor','m')

Y=X2-X1;
Compas3d(Y(1),Y(2),Y(3),'g')

%risulta
%w.'*Y=0
