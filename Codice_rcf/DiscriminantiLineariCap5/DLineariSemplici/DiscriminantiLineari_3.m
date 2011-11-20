% DiscriminantiLineari.m
%
% Discriminanti Lineari caso Wo = 0
%       t
% g(x)=W x + Wo
%
% Per un problema a 2 classi ottengo il seguente comportamento:
%
% 1) g(x) > 0 --> x è in C1 
%
% 2) g(x) < 0 --> x è in C2
%
% 3) g(x) = 0 --> x può appartenere sia a C1 che C2 : non prendo decisione
%                 IPERPIANO DI SEPARAZIONE (quando g(x) è LINEARE)
%
% Procedimento:
%
% 1) Costruisco il classificatore Lineare g(x) a partire dai seguenti
%    valori:
%
%    - w = [2 ;3;1]
%
%    - wo = 0 
%
% 2) Disegno il Piano di Separazione determinato da g(x)
%
% 3) Genero un in insieme di Pattern 3-d
%
% 4) Seleziono random 10 Pattern e li classifico
clc;clear;close all
% carattere di newline
s=sprintf('\n');

%Vettore dei Pesi
%
w=[2 ;3;1];
wo=0;
display('Vettore w (normale al piano di Separazione)')
disp(w)
display('Bias')
disp(wo)
%                            2
% Discriminante Lineare_ g: R -> R  
G=@(w,x,wo) w.'*x+wo

[x1 x2 x3]=ndgrid(-10:10);
[m n p]=size(x1);

D=[x1(:).'; x2(:).'; x3(:).'];
z=G(w,D,wo);

syms x1 x2 x3
x=[x1;x2;x3];

G=w.'*x+wo;
Sol=solve(G,'x3');
ezsurf(Sol);hold on;
N=null(w.','r')

a=N(:,1);
b=N(:,2);
compas3d(a(1),a(2),a(3),'r')
compas3d(b(1),b(2),b(3),'r')
compas3d(w(1),w(2),w(3),'m')

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
    else
        %Il punto appartiene alla Regione R2: lato Negativo
        plot3(D(1,I),D(2,I),D(3,I),'ob','MarkerSize',5,'MarkerFaceColor','b')
    end
end

%Disegno gli assi cartesiani
compas3d(10,0,0,'y')
compas3d(0,10,0,'y')
compas3d(0,0,10,'y')


%Disegno 2 punti X1 e X2 appartenenti all'Iperpiano H di separazione
A=subs(Sol,{x1,x2},{1,1});
X1=[1;1;A];
plot3(X1(1),X1(2),X1(3),'om','MarkerSize',5,'MarkerFaceColor','m')
B=subs(Sol,{x1,x2},{2,3});
X2=[2;3;B];
plot3(X2(1),X2(2),X2(3),'om','MarkerSize',5,'MarkerFaceColor','m')

Y=X2-X1;
compas3d(Y(1),Y(2),Y(3),'g')

%risulta
%w.'*Y=0
