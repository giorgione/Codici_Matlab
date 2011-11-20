% DiscriminantiLineari.m
%
% Discriminanti Lineari caso Wo > 0
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
%    - wo = 8 
%
% 2) Disegno il Piano di Separazione determinato da g(x)
%
% 4) Genero un in insieme di Pattern 3-d
%
% 5) Seleziono random 10 Pattern e li classifico
clc;clear;close all

% carattere di newline
s=sprintf('\n');

%Vettore dei Pesi
%
w=[2 ;3;1];
wo=8;
display('Vettore w (normale al piano di Separazione)')
disp(w);
display('Bias')
disp(wo)
%                            2
% Discriminante Lineare_ g: R -> R  
G=@(w,x,wo) w.'*x+wo;

[x1 x2 x3]=meshgrid(-10:10);
[m n p]=size(x1);

D=[x1(:).'; x2(:).'; x3(:).'];

% Seleziona 10 punti Random tra quelli generati con ndgrid e li classifico
% mediante il classificatore Lineare
indici=round(1+rand(1,10)*m*n*p);
TrainingData=D(:,indici);


%Valuto la funzione G(x) nei punti generati con ndgrid
z=G(w,TrainingData,wo);

syms x1 x2 x3
X=[x1;x2;x3];

F=w.'*X+wo;
Sol=solve(F,'x3');
ezsurf(Sol);

title(['Superfice di Decisione : wo=' num2str(wo)])
xlabel('x1');
ylabel('x2');
zlabel('x3');
hold on;

n=length(z);

for I=1:n
    if z(I)>0
        %Il punto appartiene alla Regione R1: lato positivo
        plot3(TrainingData(1,I),TrainingData(2,I),TrainingData(3,I),'or','MarkerSize',5,'MarkerFaceColor','r')
        %calcolo la Distanza di X dal piano
        R(I)=z(I)/norm(w);
        %P=[TrainingData(:,I) R*(w/norm(w));];
        %Disegno la proiezione di x sul piano
        %l=line(P(1,:),P(2,:),P(3,:));
        %set(l,'Color','m','LineWidth',1)
    else
        %Il punto appartiene alla Regione R2: lato Negativo
        plot3(TrainingData(1,I),TrainingData(2,I),TrainingData(3,I),'ob','MarkerSize',5,'MarkerFaceColor','b')
         %calcolo la Distanza di X dal piano
        R(I)=z(I)/norm(w);
        %P=[TrainingData(:,I) R*(w/norm(w))];
        %Disegno la proiezione di x sul piano
        %l=line(P(1,:),P(2,:),P(3,:));
        %set(l,'Color','m','LineWidth',1)
    end
end

%Disegno gli Assi
Compas3d(10,0,0,'y')
Compas3d(0,10,0,'y')
Compas3d(0,0,10,'y')

%Disegno La normale al Paino Separatore
Compas3d(w(1),w(2),w(3),'m')

%Disegno 2 punti X1 e X2 appartenenti all'Iperpiano H di separazione
A=subs(Sol,{x1,x2},{1,1});
X1=[1;1;A]
plot3(X1(1),X1(2),X1(3),'om','MarkerSize',5,'MarkerFaceColor','m')
B=subs(Sol,{x1,x2},{2,3})
X2=[2;3;B]
plot3(X2(1),X2(2),X2(3),'om','MarkerSize',5,'MarkerFaceColor','m')

Y=X2-X1;
Compas3d(Y(1),Y(2),Y(3),'g')

%risulta
%w.'*Y=0

%Calcolo la Distanza del Piano dall' Origine:
%
d=wo/norm(w);
D=(-w/norm(w))*d;

Compas3d(D(1),D(2),D(3),'r')

