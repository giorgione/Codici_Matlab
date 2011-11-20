%Derivata di un vettore numerico: Vettore delle derivate parziali
%Per calcolare le derivate di un vettore numerico di Rn basta considerare
%le sue componenti separatamente
clc;clear;

syms a1 a2 a3 a4 b1 b2;
a=[a1;a2];
b=[b1;b2];

I=eye(2);
I(:,1)*a(1)+I(:,2)*a(2)
display('Derivata del vettore a')
%Da risulta uguale ad 1
Da=diag(jacobian(a,a))

%Calolo la derivata del prodotto riga-colonna c=b'*b rispetto b e verifico che:
%
% Dc = 2*b
display('c=b.''*b')
c=b.'*b
display('Derivata di c rispetto b')
Dc=jacobian(c,b)


%Calolo la derivata del prodotto riga-colonna c=a'*b rispetto b e verifico che:
%
% Dc = a.'
display('c=a.''*b')
c=a.'*b
display('Derivata di c rispetto b')
Dc=jacobian(c,b)

%Calolo la derivata del prodotto riga-colonna c=b'*a   rispetto b e verifico che:
%
% Dc = a.'
display('c=b.''*a')
c=b.'*a
display('Derivata di c rispetto b')
Dc=jacobian(c,b)

%Quindi D(c) = D(b'*a) = D(a.'*b) =a.'

%Calolo la derivata del prodotto riga-colonna c=(b.'*a)^2   rispetto b e verifico che:
%
% Dc = 2*(b.'*a)*a.'
display('c=(b.''*a)^2')
c=(b.'*a)*(b.'*a)
display('Derivata di c rispetto b')
Dc=jacobian(c,b)
Dc1=2*(b.'*a)*a.'

%Calolo la derivata del prodotto riga-colonna c=(a.'*b)^2   rispetto b e verifico che:
%
% Dc = 2*(a.'*b)*a.'
display('c=(a.''*b)^2')
c=(a.'*b)*(a.'*b)
display('Derivata di c rispetto b')
Dc=jacobian(c,b)
Dc1=2*(a.'*b)*a.'

%Calolo la derivata del prodotto riga-colonna c=b.'A*b rispetto b con A matrice Simmetrica
%e verifico che:
%
% Dc = 2*A*b
A=[a1 0;0 a3];
c=b.'*A*b;
Dc=jacobian(c,b).';
Dc1=2*A*b;

display('Matrice A')
disp(A)
display('c=b.''*A*b')
disp(c)
display('Derivata di c rispetto b con jacobian')
disp(Dc)
display('Derivata di c rispetto b mediante formula: 2*A*b')
disp(Dc1)


%Calolo la derivata del prodotto riga-colonna c=b.'A*b rispetto b con A matrice Simmetrica
%e verifico che:
%
% Dc = 2*A*b
A=[a1 a2;a3 a4];
c=b.'*A*b;
Dc=jacobian(c,b).';
Dc1=2*A*b;

display('Matrice A')
disp(A)
display('c=b.''*A*b')
disp(c)
display('Derivata di c rispetto b con jacobian')
disp(Dc)
display('Derivata di c rispetto b mediante formula: 2*A*b')
disp(Dc1)


%Calolo la derivata del prodotto riga-colonna c=A*b rispetto b con A
%matrice generica
%e verifico che:
%
% Dc = A.'
A=[a1 a2;a3 a4];
c=A*b;
Dc=jacobian(c,b).';
Dc1=A.';

display('Matrice A')
disp(A)
display('c=A*b')
disp(c)
display('Derivata di c rispetto b con jacobian')
disp(Dc)
display('Derivata di c rispetto b mediante formula: A')
disp(Dc1)


%Calolo la derivata del prodotto riga-colonna c=A.'*b rispetto b con A
%matrice generica
%e verifico che:
%
% Dc = A
c=A.'*b;
Dc=jacobian(c,b).'
Dc1=A;

display('Matrice A')
disp(A)
display('c=A.''*b')
disp(c)
display('Derivata di c rispetto b con jacobian')
disp(Dc)
display('Derivata di c rispetto b mediante formula: A.''')
disp(Dc1)

%Calolo la derivata della norma del prodotto Matrice Vettore c=A*b rispetto b con A
%matrice generica
%e verifico che:
%
% Dc = A

c=(A*b).'*(A*b);
%applicando la proprietà di derivazione di un prodotto ho
% A.'*(A*b)+A.'*A*b = 2*A.'*A*b
Dc=jacobian(c,b).';
Dc1=2*(A.'*A)*b;
Dc=expand(Dc);
Dc1=expand(Dc1);

display('Matrice A')
disp(A)
display('c=|| (A*b) || ^2')
disp(c)
display('Derivata di c rispetto b con jacobian')
disp(Dc)
display('Derivata di c rispetto b mediante formula: 2*A.''*A*b ')
disp(Dc1)


%Calolo la derivata della norma del prodotto Matrice Vettore c=b.'*A rispetto b con A
%matrice generica
%e verifico che:
%
% Dc = A
display('c=|| (b.''*A) || ^2')
c=(b.'*A).'*(b.'*A)
%applicando la proprietà di derivazione di un prodotto ho
% A.'*(A*b)+A.'*A*b = 2*A.'*A*b
Dc=jacobian(c,b).';
Dc=expand(Dc);
Dc1=2*(A.'*A)*b;
Dc1=expand(Dc1);

display('Matrice A')
disp(A)
display('c=|| (A*b) || ^2')
disp(c)
display('Derivata di c rispetto b con jacobian')
disp(Dc)
display('Derivata di c rispetto b mediante formula: 2*A.''*A*b ')
disp(Dc1)

