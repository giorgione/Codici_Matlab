% DistribuzioneMultivariateSimboliche1_sc.m
% 
% 1) Calcolo Simbolico Distribuzione Gaussiana Multivariata per variabili 
% statisticamente Indipendenti
%
% 2) Calcolo l' integrale e Verifico che � 1
%
% 3)Calcolo E[X] e Var[X] attraverso integrazione 


clc;clear;
syms x  Ux Ox real;
X=[x];
M=[Ux];
%Matrice di Varianza Covarianza
S=[Ox^2 ];

P=GaussianaMulti(S,M,X);

%Considero la Distribuzione Normale N(3,4) <--> Ux=3 Ox=2
P=subs(P,{Ux,Ox},{3,2});
P=simplify(P);
pretty(P);

%Calcolo l' integrale e verifico che esso è 1 (sempre uguale ad 1)
A=double(int(P,x,-inf,+inf));
display('Area della Distribuzione di Probabilita N(U=3,O=2)')
disp(A);

%Calcolo la Speranza Matematica E[X]
U=double(int(x*P,x,-inf,+inf));
display('Speranza Matematica della Distribuzione di Probabilita N(U=3,O=2)')
disp(U);

%Calcolo la Varianza Var[X]=E[(X-U)^2]
Var=double(int((x-U)^2*P,x,-inf,+inf));
display('Varianza della Distribuzione di Probabilita N(U=3,O=2)')
disp(Var);


%Calcolo la probability Density function in x=[3,4]
A=double(int(P,x,3,4));
display('P( 3<= X <=4)')
disp(A);

%Calcolo la probability Density function
A=double(int(P,x,-4,4));
display('P( 4<= X <=4)')
disp(A);