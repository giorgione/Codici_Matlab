% SintesiDati_4.m
%
% Relazioni tra PCA e Autovalori ed Autovettori
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
Varianza=std(X,0,2).^2;

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
%   - S.^2: matrice degli Autovalori di  X*X.' gli autovalori sono gi�
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

display('Autovetttori della Matrice di Scatter S calcolati Mediante Svd di Xn')
disp(Vl)
disp(s)

display('Autovetttori della Matrice di Scatter S calcolati mediante eig')
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
B=Media*ones(1,N);
a=e.'*(X-B);

display('Coefficienti a(k): ')
disp(a)

% Estraggo i coefficienti(Componenti Principali delle osservazioni lungo la 
% direzione specificate dalla direzione di massimo modulo) dalla Matrice S*V della SVD
[U,S,V]=svd(Xn);
M=S*V.';
display('Coefficienti a(k): ')
disp(M(1,:))

%Per le propriet� della SVD risulta che a(k)