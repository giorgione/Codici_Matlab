% MetodoLagrange3.m
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
%  1) x^2+Y^2<=3  -> circonferenza di raggio 3 centrata nell' origine
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
plot(xt,yt,'Color','b','LineWidth',2);hold on

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
J=jacobian(L).';
F=J;
%in questo modo risolvo il sistema verificando sempre la Condizione di
%Kuhn-Tucker
F(3)=x3*G;

%Risolvo analiticamente il Sistema di Equazioni: 
Sol=solve(F(1),F(2),F(3));
Punti=[Sol.x1 Sol.x2 Sol.x3];
[m,n]=size(Punti);

disp(s)
display('Punti Candidati Estremi Locali (x,y,teta):')
disp(Punti)
h=ezplot('-x',[-3,3]);set(h,'Color','r','LineWidth',2);hold on

%Prendo solo la prima terna poichè le ultime due nn soddisfano il requisito
%di non Negatività
set(0,'CurrentFigure',fig1)
%plot(x(:,1),x(:,2),'or','MarkerSize',5,'MarkerFaceColor','g')
ezplot('-x',[-3,3]);

% Soluzione=[];
% %Vado a Verificare le condizioni di Kuhn Tucker sulle soluzioni
% for i=1:m
%     %Verifico le Condizioni di:
%     %
%     % - Non Negatività : elimino i punti Candidati che presentano un valore Negativo di Teta
%     %
%     % - Slater : G(Punti)<0 per ogni funzione Vincolo
%     Teta=double(Punti(i,3));
%     S=double(subs(G,{x1,x2},{Punti(i,1),Punti(i,2)}))
%     if( Teta>=0 && S<0 )
%         Soluzione=[Soluzione;Punti(i,:)];
%     end
% end
% display('Punti Soluzione (x,y,teta):')
% disp(Soluzione);
% disp(s);
% 
% x=Soluzione;
% if( isempty(Soluzione)==0)
%     plot3(double(x(:,1)),double(x(:,2)),double(subs(F1,{x1,x2},{x(:,1),x(:,2)})),'om','MarkerSize',5,'MarkerFaceColor','m');
% end

%Considero un punto di Minimo Soluzione
X=[-1;1;0];

%Calcolo lo jacobiano della Lagrangiana
J=jacobian(L,[x1,x2,x3]);

display('Jacobian(L(-1,1,0)):')
S=double(subs(J,{x1,x2,x3},{X(1),X(2),X(3)}));
disp(S)
disp(s)

%Calcolo la Hessiana della Lagrangiana
H=hessian(L,[x1,x2,x3]);
H=double(subs(H,{x1,x2,x3},{X(1),X(2),X(3)}));

display('Hessian(L(-1,1,0))):')
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
S=double(subs(detH,{x1,x2,x3},{X(1),X(2),X(3)}));
Val=double(subs(F1,{x1,x2},{X(1),X(2)}));
    if(S<0)
        if(H(1,2)<0 && H(1,3)<0 && H(2,3)<0)
            display(['F(' num2str(X(1)) ' , ' num2str(X(2)) ')=' num2str(Val) ' : punto di Minimo della Lagrangiana (' num2str(S) ')'])
        else
            display(['F(' num2str(X(1)) ' , ' num2str(X(2)) ')=' num2str(Val) ' : punto di Sella della Lagrangiana (' num2str(S) ')'])
        end
    end       
    if (S>0 && H(1,2)>0 && H(1,3)>0 && H(2,3)>0)
         display(['F(' num2str(X(1)) ' , ' num2str(X(2)) ')=' num2str(Val) ' : punto di Massimo della Lagrangiana  (' num2str(S) ')'])
    end
       
     







