% MixtureDensity_ML_1.m
%
% Nb: L'idea di base è che ho P(x,z| w) trattabile
%
% Data p,una Mistura di 2 Gaussiane 1-D ed un set di Dati 1-D applico l'
% algoritmo ML e la sua versione ITERATIVA per la stima del Vettore 
% delle MEDIE.
%
% p=Pw1*GaussianaMulti(1,u(1),x) + Pw2*GaussianaMulti(1,u(2),x)
%
% con:
%
% - Pw1=1/3 
%
% - Pw2=2/3
%
% - O= [1;1] -> Vettore delle Varianze Noto
%
% - u= vettore delle Medie da stimare
clc;close all;clear;

Dati=[0.608; -1.59; 0.235; 0.3949; -2.249; 2.704; -2.473; 0.672;
      0.262; 1.072; -1.773; 0.537; 3.24; 2.4; -2.499; 2.608;
      -3.458; 0.257; 2.569; 1.415; 1.410; -2.653; 1.396; 3.286; -0.712];
      

 syms x u1 u2;
 
 Pw1=1/3;
 Pw2=2/3;
 
 %Costruisco la Mistura
 Fun_Mixture=@(x,u) Pw1*GaussianaMulti(1,u(1),x) + Pw2*GaussianaMulti(1,u(2),x);
 
 n=length(Dati);
 
 L=0;
 %Calcolo L la Log-Likelyhood function
 for i=1:n
     L=L+log(Fun_Mixture(Dati(i),[u1;u2]));
 end
 
 ezsurf(L,[-20,20,-20,20])
 f1 =gcf;
 %Definisco i settaggi dell' OTTIMIZAZIONE
 options = optimset('Display','final','MaxFunEvals',1500,'MaxIter',1500,'TolFun',1e-10,'TolX',1e-10);
 
 f2=figure;
 ezcontour(-L,[-10,10,-10,10],120);hold on;
 
 Uo=ginput(1);
 plot(Uo(1),Uo(2),'or','MarkerFaceColor','r')
 
 %Calcolo il Massimo della Funzione Verosimiglianza per la Mistura
 fun=@(u) -Mixture_LogLike([1;1],u,Dati,Pw1,Pw2);
 [Media1 ,val,exitflag,output]=fminsearch(fun,Uo,options);
 
 plot(Media1(1),Media1(2),'ob','MarkerFaceColor','b')
 
 set(0,'CurrentFigure',f1);hold on
 plot3(Uo(1),Uo(2),-fun(Uo),'or','MarkerFaceColor','r')
 plot3(Media1(1),Media1(2),-val,'ob','MarkerFaceColor','b');
 
 set(0,'CurrentFigure',f2);hold on
 
 display('Metodo Iterativo: premere un tasto per continuare..')
 pause;
 %Applico l' algoritmo Iterativo per la Stima del Vettore delle Medie:
 %
 U=Uo.';
 plot(Uo(1),Uo(2),'or','MarkerSize',2,'MarkerFaceColor','r')
 for J=1:30
     U(:,J+1)=0;
     Den=0;
     P1=zeros(1,25);
     P2=zeros(1,25);
     
     %% E-Step:  
     %Per ogni singolo campione vado a Misurare le Responsabilities
     for K=1:25 
         
        % Valuto la Mistura nel Punto x(k)  P(x(k) | �)
        p(K)=Mixture([1;1],U(:,J),Dati(K),Pw1,Pw2);
        
        
        %Calcolo la Verosimiglianza del Campione Xk alla componente 1 della
        %Mistura: P(w(i) | x(k) , �(i))
        %
        % P(w(1) | x(k) , �(1))
        p1(K)=GaussianaMulti(1,U(1,J),Dati(K));
        
        %Calcolo la Verosimiglianza del Campione Xk alla componente 2 della
        %Mistura: P(w(i) | x(k) , �(i))
        %
        % P(w(2) | x(k) , �(2))
        p2(K)=GaussianaMulti(1,U(2,J),Dati(K));
        
        %Responsabilities: misura quanto ciascuna componente explains il
        %dato osservato
        P1(K)=p1(K)*Pw1/p(K);
        P2(K)=p2(K)*Pw2/p(K);
        
     end
     %% M-Step
     %Aggiorno le Medie
     U(1,J+1)=(P1*Dati)/sum(P1);
     U(2,J+1)=(P2*Dati)/sum(P2);
     
     
     
     plot(U(1,J+1),U(2,J+1),'om','MarkerSize',2,'MarkerFaceColor','m')
     %Disegna la linea tra il nuovo ed il Vecchio Punto
     l=line( [U(1,J+1);U(1,J)] , [U(2,J+1);U(2,J)]);
     set(l,'Color','m','LineWidth',2);
     pause;
 end
 plot(U(1,end),U(2,end),'og','MarkerSize',2,'MarkerFaceColor','g');
 
 set(0,'CurrentFigure',f1)
 plot3(U(1,end),U(2,end),subs(L,{u1,u2},{U(1,end),U(1,end)}),'og','MarkerSize',5,'MarkerFaceColor','g');
 
 Media2=U(:,end);

 figure;
 
 syms x
 F=Fun_Mixture(x,Media1);
 ezplot(F);
 hold on;
 plot(Dati,zeros(1,25),'og','MarkerSize',2,'MarkerFaceColor','g');
 plot(Media2(1),0,'or','MarkerSize',2,'MarkerFaceColor','r')
 plot(Media2(2),0,'or','MarkerSize',2,'MarkerFaceColor','r')
 title('Mistura di Gaussiane con Media Stimata mediante ML')