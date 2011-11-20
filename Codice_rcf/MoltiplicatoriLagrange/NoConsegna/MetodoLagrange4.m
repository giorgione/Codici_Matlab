% MetodoLagrange4.m
%
%Risolvo attraverso il metodo di Lagrange il seguente problema di
%ottimizazione a vincoli non lineari (DISUGUAGLIANZE) in R^2
%
% Funzione obbiettivo:
%
%  F(x,y)= cos(x1)*sin(x2)
%
% Vincoli non Lineari (RICERCA SUI PUNTI DI FRONTIERA ed INTERNI <--> G(X)<=0)
%
%
%  1) x^2+y^2-9<=0  -> circonferenza di raggio 3 centrata nell' origine
%
%  2) (x1-pi)^2+(x2-(-pi))^2-r2^2 <=0 

clc;clear;close all
% carattere di newline
s=sprintf('\n');

syms x1 x2 x3 x4 real
r1=4;
r2=3;
F1=x1^2+x2^2+2*x1*x2;
G1=x1^2+x2^2-r1^2;
G2=(x1-pi)^2+(x2-(-pi))^2-r2^2;

%Disegno la curva definita dai vincoli sul grafico di F
[xt1,yt1]=DisegnaCirconferenza([0;0],r1,30,0);
fig1=figure;
plot(xt1,yt1,'Color','b','LineWidth',2);hold on;

[xt2,yt2]=DisegnaCirconferenza([pi;-pi],r2,30,0);
plot(xt2,yt2,'Color','b','LineWidth',2);

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

%Risolvo il Problema di Lagrange per il Calcolo dei Punti Stazionari della
%Lagrangiana  -> sistema di Equazioni Non Lineari risolte mediante il
%metodo di Newton con punto iniziale Xo
J=jacobian(L).';
F=J;
%in questo modo risolvo il sistema verificando sempre la Condizione di
%Kuhn-Tucker
F(3)=x3*G1;
F(4)=x4*G2;


%Risolvo analiticamente il Sistema di Equazioni: 
Sol=solve(F(1),F(2),F(3),F(4),'x1','x2','x3','x4');
Punti=[Sol.x1 Sol.x2 Sol.x3 Sol.x4];
[m,n]=size(Punti);

disp(s)
display('Punti Candidati Estremi Locali (x,y,teta):')
disp(Punti)

Soluzione=[];
%Vado a Verificare le condizioni di Kuhn Tucker sulle soluzioni
for i=1:m
    %Verifico le Condizioni di:
    %
    % - Non Negatività : elimino i punti Candidati che presentano un valore Negativo di Teta
    %
    % - Slater : G(Punti)<0 per ogni funzione Vincolo
    Teta1=double(Punti(i,3));
    Teta2=double(Punti(i,4));
    %S1=double(subs(G1,{x1,x2},{Punti(i,1),Punti(i,2)}));
    %S2=double(subs(G2,{x1,x2},{Punti(i,1),Punti(i,2)}));

    if( Teta1>=0 && Teta2>=0)
        Soluzione=[Soluzione;Punti(i,:)];
    end
end
display('Punti Soluzione (x,y,teta):')
disp(Soluzione);
disp(s);

h=ezplot('-x',[-3,3]);set(h,'Color','r','LineWidth',2);hold on

%Prendo solo la prima terna poichè le ultime due nn soddisfano il requisito
%di non Negatività
set(0,'CurrentFigure',fig1)
%plot(x(:,1),x(:,2),'or','MarkerSize',5,'MarkerFaceColor','g')
ezplot('-x',[-3,3]);

% x=Soluzione;
% if( isempty(Soluzione)==0)
%     plot3(double(x(:,1)),double(x(:,2)),double(subs(F1,{x1,x2},{x(:,1),x(:,2)})),'om','MarkerSize',5,'MarkerFaceColor','m');
%     set(0,'CurrentFigure',fig1)
%     plot(x(:,1),x(:,2),'or','MarkerSize',5,'MarkerFaceColor','g')
% end


%Considero un punto di Minimo Soluzione
X=[-1;1;0;0];

%Calcolo lo jacobiano della Lagrangiana
J=jacobian(L,[x1,x2,x3,x4]);

display('Jacobian(L(-1,1,0,0)):')
S=double(subs(J,{x1,x2,x3,x4},{X(1),X(2),X(3),X(4)}));
disp(S)
disp(s)

%Calcolo la Hessiana della Lagrangiana
H=hessian(L,[x1,x2,x3,x4]);
H=double(subs(H,{x1,x2,x3,x4},{X(1),X(2),X(3),X(4)}));

display('Hessian(L(-1,1,0,0))):')
disp(H)
disp(s)

%Condizioni da Verificare per il punto di Massimo-Minimo-Sella:
%
% Sia H la matrice Hessiana
%
% Xo ->MASSIMO se :
%
%       - det(H) > 0
%
%       - H(i,j) > 0  per ogni i!=j
%
% Xo ->MINIMO se :
%
%       - det(H) < 0
%
%       - H(i,j) < 0  per ogni i!=j
%
% Xo -> SELLA se:
%
%       - det(H) < 0
%Calcolo il Determinante
detH=det(H);
S=double(subs(detH,{x1,x2,x3,x4},{X(1),X(2),X(3),X(4)}));
Val=double(subs(F1,{x1,x2},{X(1),X(2)}));
    if(S<0)
        if(H(1,2)<0 && H(1,3)<0 && H(1,4)<0 && H(2,3)<0 && H(2,4)<0 && H(3,4)<0)
            display(['F(' num2str(X(1)) ' , ' num2str(X(2)) ')=' num2str(Val) ' : punto di Minimo della Lagrangiana (' num2str(S) ')'])
        else
            display(['F(' num2str(X(1)) ' , ' num2str(X(2)) ')=' num2str(Val) ' : punto di Sella della Lagrangiana (' num2str(S) ')'])
        end
    end       
    
    if (S>0 && H(1,2)>0 && H(1,3)>0 && H(1,4)>0 && H(2,3)>0 && H(2,4)>0 && H(3,4)>0)
         display(['F(' num2str(X(1)) ' , ' num2str(X(2)) ')=' num2str(Val) ' : punto di Massimo della Lagrangiana  (' num2str(S) ')'])
    end

