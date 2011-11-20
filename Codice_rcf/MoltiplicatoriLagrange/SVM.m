clc;clear;close all;
k=menu('Seleziona i Dati:','1)Dati L.S. mediante retta non per l'' origine','2) Dati L.S. mediante retta per l''origine');
switch k
        case 1
        % Il Metodo del Percettrone BIASED  converge alla Soluzione 
        % ottima poichè esiste una retta di Separazione passante per
        % l' origine che produce 0 errori di classificazione sul TRAINING
        % SET.
        x1=5+randn(1,10);
        y1=5+randn(1,10);

        x2=randn(1,10);
        y2=randn(1,10);
        
        %Insericsco i dati in un unico Vettore
        Vet=[x1 x2;y1 y2];
        n=length(x1);
        
    case 2
        %Il Metodo del Percettrone NOT BIASED converge alla soluzione OTTIMA
        % con 0 errori di classificazione sul TRAINING SET.
        x1=-2+randn(1,10);
        y1=-2+randn(1,10);

        x2=2+randn(1,10);
        y2=2+randn(1,10);
        
        %Insericsco i dati in un unico Vettore
        Vet=[x1 x2;y1 y2];
        n=length(x1);
end

%Suppongo che:
%- i primi 5 vettori (y1...y5) appartengono ad W1 --> label 1
%- gli altri 5 (y6..,y10) appartengano ad W2 --> label -1
x=Vet;
%V=[ones(1,2*n);V]
%V(:,n+1:2*n)=-V(:,n+1:2*n);

plot(Vet(1,1:n), Vet(2,1:n),'or','MarkerSize',5,'MarkerFaceColor','r');
hold on;
plot(Vet(1,n+1:2*n), Vet(2,n+1:2*n),'ob','MarkerSize',5,'MarkerFaceColor','b')
plot(0,0,'dg','MarkerSize',3,'MarkerFaceColor','g')

%Numero di Campioni
N=2*n;
y=ones(N,1);
y(1:N/2)=-1;
B=zeros(N,N);

for i=1:N
    for j=1:N
        B(i,j)=(y(i)*y(j))*x(:,i).'*x(:,j);
    end
end

options = optimset('Display','iter','LargeScale','off','TolCon',1e-16);
L=@(x) Lagrangiana_SVM(x,B);
V_L=@(x)Vincoli_Lagrangiana(x,y);
%[a,fval] = fmincon(L,rand(N,1),[],[],[],[],[],[],V_L,options);


f=ones(N,1);
A=-eye(N);
b=zeros(N,1);
Aeq=y.';
beq=0;
[a,fval]=quadprog(B,-f,A,b,Aeq,beq,0,20,[],options);
pause;
%Calcolo il vettore w
w=x*(a.*y);

%Calcolo il bias b
for i=1:N
    b(i)=a(i)-a(i)*y(i)*x(:,i).'*w;
end
b=mean(y.'*b);

syms X Y real;
S=solve([b;w].'*[1;X;Y],Y);
ezplot(S,[-10 10]);axis tight
title('Retta di Separazione a Massimo Margine');



