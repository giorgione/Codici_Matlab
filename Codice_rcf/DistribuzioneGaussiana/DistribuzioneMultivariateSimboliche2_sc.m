%Calcolo Simbolico Distribuzione Gaussiana bivariata 
clc;clear;
syms x y  Ux Uy Ox Oy Oxy real;
X=[x;y];
M=[Ux;Uy];
%Matrice di Varianza Covarianza
S=[Ox^2 Oxy;
   Oxy  Oy^2];

P=GaussianaMulti(S,M,X);

%Considero la Distribuzione Normale N(3,2)
P=subs(P,{Ux,Uy,Ox,Oy,Oxy},{3,1,1,2,0});
P=simplify(P);
pretty(P);

%Calcolo l' integrale e verifico che esso � 1 (sempre uguale ad 1)
I1=int(P,x,-inf,+inf);
A=double(int(I1,y,-inf,+inf));
display('Area della Distribuzione di Probabilita N(U=[3;1],O=[1 0;0 4])')
disp(A);

%Calcolo la Speranza Matematica E[X] come Integrale doppio risolto prima
%rispetto ad x e poi rispetto a y
U1=int(X*P,x,-inf,+inf);
U=double(int(U1,y,-inf,+inf));
display('Speranza Matematica della Distribuzione di Probabilita N(U=[3;1],O=[1 0;0 4])')
disp(U);

%Calcolo la Speranza Matematica E[X] come Integrale doppio risolto prima
%rispetto ad y e poi rispetto a x
U1=int(X*P,y,-inf,+inf);
U=double(int(U1,x,-inf,+inf));
display('Speranza Matematica della Distribuzione di Probabilita N(U=[3;1],O=[1 0;0 4])')
disp(U);


%Calcolo la Matrice di Varianza Covarianza Var[X]=E[(X-U)^2] come Integrale 
%doppio risolto prima rispetto ad x e poi rispetto a y
Var1=int(((X-U)*(X-U).')*P,x,-inf,+inf);
Var=double(int(Var1,y,-inf,+inf));

ezsurf(P)
display('Matrice di Varianza-Covarianza della Distribuzione di Probabilita N(U=[3;1],O=[1 0;0 4])')
disp(Var);

%Calcolo la Matrice di Varianza Covarianza Var[X]=E[(X-U)^2] come Integrale 
%doppio risolto prima rispetto ad y e poi rispetto a x
Var1=int(((X-U)*(X-U).')*P,y,-inf,+inf);
Var=double(int(Var1,x,-inf,+inf));
display('Matrice di Varianza-Covarianza della Distribuzione di Probabilita N(U=[3;1],O=[1 0;0 4])')
disp(Var);