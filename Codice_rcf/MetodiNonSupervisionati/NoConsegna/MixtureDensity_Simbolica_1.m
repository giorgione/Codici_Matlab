% MixtureDensity_Simbolica.m
% Costruisco Simbolicamente la Mistura di Densità
clc;close all;clear;

Dati=[0.608; -1.59; 0.235; 0.3949; -2.249; 2.704; -2.473; 0.672;
      0.262; 1.072; -1.773; 0.537; 3.24; 2.4; -2.499; 2.608;
      -3.458; 0.257; 2.569; 1.415; 1.410; -2.653; 1.396; 3.286; -0.712];
      

 syms x u1 u2 o1 o2 Pw1 Pw2;
 
 %Pw1=1/3;
 %Pw2=2/3;
 
 %Costruisco la Mistura
 Fun_Mixture=@(x,u,o,P) P(1)*GaussianaMulti(o(1),u(1),x) + P(2)*GaussianaMulti(o(2),u(2),x);
 
 n=length(Dati);
 
 L=0;
 %Calcolo L la Log-Likelyhood function
 for i=1:n
     L=L+log(Fun_Mixture(Dati(i),[u1;u2],[o1;o2],[Pw1,Pw2]));
 end
 
 %ezsurf(L,[-20,20,-20,20])
 %f1 =gcf;
 %Definisco i settaggi dell' OTTIMIZAZIONE
 options = optimset('Display','final','MaxFunEvals',1500,'MaxIter',1500,'TolFun',1e-10,'TolX',1e-10);
 
 %f2=figure;
 %ezcontour(-L,[-10,10,-10,10],120);hold on;
 
 %Uo=ginput(1);
 %plot(Uo(1),Uo(2),'or','MarkerFaceColor','r')
 Xo=randn(1,4);
 Xo=[Xo .5 .5];
 
 %Calcolo il Massimo della Funzione Verosimiglianza per la Mistura
 fun=@(X) -Mixture_LogLike(X(1:2),X(3:4),Dati,X(5),X(6));
 [x ,val,exitflag,output]=fminsearch(fun,Xo,options)
 
%  plot(x(1),x(2),'ob','MarkerFaceColor','b')
%  
%  set(0,'CurrentFigure',f1);hold on
%  plot3(Uo(1),Uo(2),-fun(Uo),'or','MarkerFaceColor','r')
%  plot3(x(1),x(2),-val,'ob','MarkerFaceColor','b');
%  
%  set(0,'CurrentFigure',f2);hold on
%  
%  display('Metodo Iterativo: premere un tasto per continuare..')
%  pause;
%  %Applico l' algoritmo Iterativo per la Stima del Vettore delle Medie:
%  %
%  U=Uo.';
%  plot(Uo(1),Uo(2),'or','MarkerSize',2,'MarkerFaceColor','r')
U=randn(2,1);
O=randn(2,1);
P=randn(2,1);
 for J=1:5
     U(:,J+1)=0;
     O(:,J+1)=0;
     P(:,J+1)=0;
     Den=0;
     P1=zeros(1,25);
     P2=zeros(1,25);
     for K=1:25 
        Den=Den+Mixture(O(:,J),U(:,J),Dati(K),P(1,J),P(2,J));
        
        %Calcolo la Verosimiglianza del Campione Xk alla componente 1 della
        %Mistura
        P1(K)=(GaussianaMulti(O(1,J).^2,U(1,J),Dati(K))*P(1,J))/Den;
        
        %Calcolo la Verosimiglianza del Campione Xk alla componente 2 della
        %Mistura
        P2(K)=(GaussianaMulti(O(2,J).^2,U(2,J),Dati(K))*P(2,J))/Den;
        
        
     end
     
        %Stimo la probabilità a priori
        P(1,J+1)=mean(P1)
        P(2,J+1)=mean(P2);
        
        %Stimo la media
        U(1,J+1)=(P1.'*Dati)/Den;
        U(2,J+1)=(P2.'*Dati)/Den;
        
        O(1,J+1)=O(1,J+1)+P1(K)*(Dati(K)- U(1,J+1))^2;
        O(2,J+1)=O(2,J+1)+P2(K)*(Dati(K)- U(2,J+1))^2;
     
     U(1,J+1)=U(1,J+1)/sum(P1);
     U(2,J+1)=U(2,J+1)/sum(P2);
     
     O(1,J+1)=O(1,J+1)/sum(P1);
     O(2,J+1)=O(2,J+1)/sum(P2);
     
    
     pause;
     
     %plot(U(1,J+1),U(2,J+1),'om','MarkerSize',2,'MarkerFaceColor','m')
     %Disegna la linea tra il nuovo ed il Vecchio Punto
     %l=line( [U(1,J+1);U(1,J)] , [U(2,J+1);U(2,J)]);
     %set(l,'Color','m','LineWidth',2);
     %pause;
 end
%  plot(U(1,end),U(2,end),'og','MarkerSize',2,'MarkerFaceColor','g');
%  
%  set(0,'CurrentFigure',f1)
%  plot3(U(1,end),U(2,end),subs(L,{u1,u2},{U(1,end),U(1,end)}),'og','MarkerSize',5,'MarkerFaceColor','g');
%  