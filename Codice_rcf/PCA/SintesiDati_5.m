% SintesiDati_5.m
%
clc;clear ;close all
% carattere di newline
s=sprintf('\n');

%Genero un un Set di N dati 3-Dimensionali 
N=10;

a=10;
b=25;
x1=a+(b-a)*rand(1,N);
x2=a+(b-a)*rand(1,N);
x3=a+(b-a)*rand(1,N);

X=[x1;x2;x3];
%Calcolo la media dei Dati: essa rappresenta una Rappresentazione Sintetica
%                           dei dati zero-dimensionale
  
Media=mean(X,2);
Varianza1=std(X,0,2).^2;
Varianza2=std(X,1,2).^2;

%Disegno i dati
plot3(X(1,:),X(2,:),X(3,:),'or','MarkerFaceColor','r');hold on;

%Disegno la Media, intesa come rappresentazione sintetica
plot3(Media(1),Media(2),Media(3),'og','MarkerFaceColor','g');
grid on;

%Calcolo la Matrice di Scatter
Scatter=MatriceScatter(X,Media)

%Matrice di Varianza-Covarianza dei Dati
B=X*X.';

%Matrice dei Dati Normalizzati: essa concide con la Matrice di Scatter
Xn=(X-Media*ones(1,N));


%Matrice di Varianza Covarianza dei dati normalizzati: risulta uguale alla
%matrice di Scatter
Bn=Xn*Xn.';

% Calcolo la svd dei dati Normalizzati per ottenere Autovalori ed Autovettori
% della Matrice d Scatter: 
%
%   - U : matrice degli Autovettori di X*X.'
%
%   - S.^2: matrice degli Autovalori di  X*X.' gli autovalori sono già
%           oridinati
%
%
[U,S,V]=svd(Xn); 
%U: autovettori di Xn*Xn.' , ortogonali e base per R(Xn) spazio delle
%osservazioni
Vl=U;
%autovalori
lamda=S.^2;


% Calcolo Autovalori ed Autovettori della Matrice d Scatter: occorre
% ordinate autovalori ed autovettori
%
[Xl l]=eig(Scatter);
%ordina gli autovalori (ed i rispettivi autovettori) in senso decrescente
d=diag(l);
[lamdaNew ind]=sort(d,'descend');
l=diag(lamdaNew);
Xl=Xl(:,ind);

display('Autovettori della Matrice di Scatter S calcolati Mediante Svd di Xn')
disp(Vl)
disp(s)

display('Autovettori della Matrice di Scatter S calcolati mediante eig')
disp(Xl)

%Calcolo la retta dei minimi quadrati e la disegno
e=Xl(:,1);
t=-10:.25:10;
n=length(t)
Retta=[];
for i=1:n
    %Equazione parametrica della Retta
    Retta=[Retta Media+t(i)*e;];
end
%Disegno la Retta, intesa come rappresentazione sintetica
plot3(Retta(1,:),Retta(2,:),Retta(3,:),'LineWidth',2);


%Calcolo i coefficienti a(k) associati alle osservazioni:
%
% a(k)= e.'*(x(k)-Media)

a=e.'*(Xn);

display('Coefficienti a(k):')
disp(a)
disp(s)

%Per le proprietà della SVD risulta che a(k)
%
%
Vettore3D_Applicato(Media,Xl(1,1),Xl(2,1),Xl(3,1),[0 1 0])
Vettore3D_Applicato(Media,Xl(1,2),Xl(2,2),Xl(3,2),[0 1 0])
Vettore3D_Applicato(Media,Xl(1,3),Xl(2,3),Xl(3,3),[0 1 0])

Vt=V.';
%Calcolo la prima osservazione x
display('Vettore x1 come combinazione lineare degli autovettori di Scatter:')
disp(U*S*Vt(:,1))
disp(s)

display('Vettore x1 originale:')
disp(Xn(:,1))
display(s)

display('Coefficienti a(k) estratti da SVD:')
disp(S(1,1)*Vt(:,1))
disp(s)

figure;
%Disegno i dati Normalizzati: Media Nulla , Varianza V
plot3(Xn(1,:),Xn(2,:),Xn(3,:),'or','MarkerFaceColor','r');hold on;
%Disegno la Media, intesa come rappresentazione sintetica
plot3(0,0,0,'og','MarkerFaceColor','g');

%Disegno le direzioni e di distanza minima
Xl=200*Xl;
Vettore3D(Xl(1,1),Xl(2,1),Xl(3,1),[0 1 1])
Vettore3D(Xl(1,2),Xl(2,2),Xl(3,2),[0 1 1])
Vettore3D(Xl(1,3),Xl(2,3),Xl(3,3),[0 1 1])
hold on
colormap;
[x y z]=ellipsoid(0,0,0,lamda(1,1),lamda(2,2),lamda(3,3));
surf(x,y,z)


