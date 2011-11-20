% MetodoLagrange2.m
%
%Risolvo attraverso il metodo di Lagrange il seguente problema di
%ottimizazione a vincoli non lineari in R^2
%
% Funzione obbiettivo:
%
%  F(x,y)= cos(x1)*sin(x2)
%
% Vincoli non Lineari (RICERCA SUI PUNTI DI FRONTIERA <--> G(X)=0):
%
%  1) x^2+y^2-9=0  -> circonferenza di raggio 3 centrata nell' origine
%
%  2) (x-2)^2+(y-1)^2-4=0 

clc;clear;close all
% carattere di newline
s=sprintf('\n');

syms x1 x2 x3 x4
r1=3;
r2=2;

F1=x1^2+x2^2+2*x1*x2;
G1=x1^2+x2^2-r1^2;
G2=(x1-2)^2+(x2-1)^2-r2^2;

%Disegno la curva definita dai vincoli sul grafico di F
[xt1,yt1]=DisegnaCirconferenza([0;0],r1,30,0);
fig1=figure;
plot(xt1,yt1,'Color','b','LineWidth',2);hold on;

[xt2,yt2]=DisegnaCirconferenza([2;1],r2,30,0);
plot(xt2,yt2,'Color','b','LineWidth',2);

Xo=ginput(4);
hold on
plot(Xo(:,1),Xo(:,2),'og','MarkerSize',5,'MarkerFaceColor','g');
Xo=[Xo randn(4,1) randn(4,1)];

xt=[xt1;xt2];
yt=[yt1;yt2];

Fnum=xt1.^2+yt1.^2+2*xt1.*yt1;
figure;
h=plot3(xt1,yt1,Fnum);
set(h,'Color','b','LineWidth',2);hold on

Fnum=xt2.^2+yt2.^2+2*xt2.*yt2;
h=plot3(xt2,yt2,Fnum);
set(h,'Color','b','LineWidth',2);hold on
ezsurf(F1,[-5 5 -5 5]);


teta=[x3;x4];
%Equazione di Lagrange con 2 moltiplicatori
L=F1+teta(1).'*G1+teta(2).'*G2;
%Vettore delle variabili:
%[  x      ]   ->  x1
%[  y      ]   ->  x2
%[  teta1  ]   ->  x3
%[  teta2  ]   ->  x4

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
    display(['Punto Iniziale : X(' num2str(Xo(I,1)) ',' num2str(Xo(I,2)) ',' num2str(Xo(I,3)) ',' num2str(Xo(I,4)) ')'])
    Nk1=0;
    Nk=1;
    k=1;
    tol=1e-4;
    N=1;
    while k<=400 
        if(N>tol)
        
            if(mod(k,40)==0)
                display('.')
            end
            
        %Valuta lo Iacobiano nel punto Xk
        J_k=subs(JF,{x1,x2,x3,x4},{X(k,1), X(k,2),X(k,3),X(k,4)});

        %Valuta la Funzione nel punto Xk
        Fun(:,k)=subs(F,{x1,x2,x3,x4},{X(k,1), X(k,2),X(k,3),X(k,4)});

        %Risolvo il sistema 
        Ck=J_k\(-Fun(:,k));

        %Applico l' iterazione del metodo di Newton
        X(k+1,:)=X(k,:)+Ck.';

        %Valuto la Funzione nel nuovo punto approssimante la Soluzione
        Fun(:,k+1)=subs(F,{x1,x2,x3,x4},{X(k,1), X(k,2),X(k,3),X(k,4)});


        %Criterio di Arresto:
        %  
        % || X(k+1,:) - X(k,:) || < tol
        %N=norm(X(k+1,:)-X(k,:),2);
        
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
    x(I,4)=double(X(end,4));
    
    display(['Estremo Locali (x,y,teta): X(' num2str(x(I,1)) ',' num2str(x(I,2)) ',' num2str(x(I,3)) ',' num2str(x(I,4)) ')'])
    disp(s)
end
disp(s)
display('Punti Estremi Locali (x,y,teta1,teta2):')
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

%Calcolo la Hessiana della Lagrangiana
H=hessian(L,[x1,x2,x3,x4])

%Calcolo il Determinante
detH=det(H);

[m,n]=size(x)

for k=1:m
    S=double(subs(detH,{x1,x2,x3,x4},{x(k,1), x(k,2),x(k,3),x(k,4)}));
    Val=double(subs(F1,{x1,x2},{x(k,1), x(k,2)}));
    if(S<0)
         display(['F(' num2str(x(k,1)) ' , ' num2str(x(k,2)) ')=' num2str(Val) ' : punto di Minimo Locale Vincolato (' num2str(S) ')'])
    else
         display(['F(' num2str(x(k,1)) ' , ' num2str(x(k,2)) ')=' num2str(Val) ' : punto di Massimo Locale Vincolato (' num2str(S) ')'])
    end 
end