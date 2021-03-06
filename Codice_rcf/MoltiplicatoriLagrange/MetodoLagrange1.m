% MetodoLagrange1.m
%
%Risolvo attraverso il metodo di Lagrange il seguente problema di
%ottimizazione a vincoli non lineari in R^2
%
% Funzione obbiettivo:
%
%  F(x,y)= x1^2+x2^2+2*x1*x2
%
% Vincoli non Lineari (RICERCA SUI PUNTI DI FRONTIERA <--> G(X)=0):
%
%  1) x^2+Y^2-9=0  -> circonferenza di raggio 3 centrata nell' origine
%

clc;clear;close all
% carattere di newline
s=sprintf('\n');

syms x1 x2 x3 
r=3;

F1=x1^2+x2^2+2*x1*x2;
G=x1^2+x2^2-r^2;

%Disegno la curva definita dai vincoli sul grafico di F
[xt,yt]=DisegnaCirconferenza([0;0],r,30,0);
fig1=figure;
plot(xt,yt,'Color','b','LineWidth',2);
Xo=ginput(4);
hold on
plot(Xo(:,1),Xo(:,2),'og','MarkerSize',5,'MarkerFaceColor','g');
Xo=[Xo randn(4,1)];

%l=ezplot(G(1));set(l,'Color','r','LineWidth',2);hold on
Fnum=xt.^2+yt.^2+2*xt.*yt;
figure;
h=plot3(xt,yt,Fnum);
set(h,'Color','b','LineWidth',2);hold on
ezsurf(F1,[-5 5 -5 5]);


teta=[x3];
%Equazione di Lagrange
L=F1+teta.'*G;
%Vettore delle variabili:
%[  x      ]   ->  x1
%[  y      ]   ->  x2
%[  teta1  ]   ->  x3

%Risolvo il Problema di Lagrange per il Calcolo dei Punti Stazionari della
%Lagrangiana  -> sistema di Equazioni Non Lineari risolte mediante il
%metodo di Newton con punto iniziale Xo
F=jacobian(L).';
JF=jacobian(F);
display('Metodo di Newton per la risoluzione del Sistema di Equazioni Non Lineari')
disp(s);

for I=1:4
    %Seleziona il punto Iniziale
    X=Xo(I,:);
    display(['Punto Iniziale : X(' num2str(Xo(I,1)) ',' num2str(Xo(I,2)) ',' num2str(Xo(I,3)) ')'])
    Nk1=0;
    Nk=1;
    k=1;
    tol=1e-16;
    N=1;
    while k<=400 
        if(N>tol)
        
            if(mod(k,40)==0)
                display('.')
            end
            
        %Valuta lo Iacobiano nel punto Xk
        J_k=subs(JF,{x1,x2,x3},{X(k,1), X(k,2),X(k,3)});

        %Valuta la Funzione nel punto Xk
        Fun(:,k)=subs(F,{x1,x2,x3},{X(k,1), X(k,2),X(k,3)});

        %Risolvo il sistema 
        Ck=J_k\(-Fun(:,k));

        %Applico l' iterazione del metodo di Newton
        X(k+1,:)=X(k,:)+Ck.';

        %Valuto la Funzione nel nuovo punto approssimante la Soluzione
        Fun(:,k+1)=subs(F,{x1,x2,x3},{X(k,1), X(k,2),X(k,3)});


        
        %Criterio di Arresto:
        %  
        % || F(k+1,:)|| < tol
        N=norm(Fun(:,k+1),2);
        
        %incrementa l' iterazione
        k=k+1;
        else
            display(['Arresto del metodo dopo ' num2str(k) ' passi N=' num2str(N)])
            break
        end 

    end
    plot3(double(X(end,1)),double(X(end,2)),double(subs(F1,{x1,x2},{X(end,1),X(end,2)})),'og','MarkerSize',5,'MarkerFaceColor','g');
    x(I,1)=double(X(end,1));
    x(I,2)=double(X(end,2));
    x(I,3)=double(X(end,3));
    display(['Estremo Locali (x,y,teta): X(' num2str(x(I,1)) ',' num2str(x(I,2)) ',' num2str(x(I,3)) ')'])
    disp(s)
end
disp(s)
display('Punti Estremi Locali (x,y,teta):')
disp(x)

set(0,'CurrentFigure',fig1)
plot(x(:,1),x(:,2),'or','MarkerSize',5,'MarkerFaceColor','g')
Punto=[Xo(1,1:2) ;x(1,1:2)];
l=line(Punto(:,1),Punto(:,2));set(l,'Color','m','LineWidth',2)

Punto=[Xo(2,1:2) ;x(2,1:2)];
l=line(Punto(:,1),Punto(:,2));set(l,'Color','m','LineWidth',2)

Punto=[Xo(3,1:2) ;x(3,1:2)];
l=line(Punto(:,1),Punto(:,2));set(l,'Color','m','LineWidth',2)

Punto=[Xo(4,1:2) ;x(4,1:2)];
l=line(Punto(:,1),Punto(:,2));set(l,'Color','m','LineWidth',2)

%Calcolo l' Hessiana della Lagrangiana
H=hessian(L,[x1,x2,x3])

%Calcolo il Determinante
detH=det(H);

[m,n]=size(x)

for k=1:m
    S=double(subs(detH,{x1,x2,x3},{x(k,1), x(k,2),x(k,3)}));
    Val=double(subs(F1,{x1,x2},{x(k,1), x(k,2)}));
    if(S<0)
         display(['F(' num2str(x(k,1)) ' , ' num2str(x(k,2)) ')=' num2str(Val) ' : punto di Minimo Locale Vincolato (' num2str(S) ')'])
    else
         display(['F(' num2str(x(k,1)) ' , ' num2str(x(k,2)) ')=' num2str(Val) ' : punto di Massimo Locale Vincolato (' num2str(S) ')'])
    end 
end