% MixtureDensity_ML_3.m
% 
%
% Data p,una Mistura di 2 Gaussiane 1-D ed un set di Dati 1-D applico l'
% algoritmo ML e la sua versione ITERATIVA per la stima del Vettore 
% delle MEDIE , delle VARIANZE e delle probabilit� a PRIORI .
%
% p=Pw1*GaussianaMulti(o(1),u(1),x) + Pw2*GaussianaMulti(o(2),u(2),x)
%
% con:
%
% - Pw1= probabilit� a PRIORI della classe w1 
%
% - Pw2= probabilit� a PRIORI della classe w2
%
% - O= [1;1] -> Vettore delle Varianze da Stimare
%
% - u= vettore delle Medie da stimare
%
% E' un problema di ottimizazione VINCOLATA dai VINCOLI:
%
% - Pw1+Pw2=1
%
% - Pw1 >=0
%
% - Pw2 >=0

clc;close all;clear;

Dati=[0.608; -1.59; 0.235; 0.3949; -2.249; 2.704; -2.473; 0.672;
      0.262; 1.072; -1.773; 0.537; 3.24; 2.4; -2.499; 2.608;
      -3.458; 0.257; 2.569; 1.415; 1.410; -2.653; 1.396; 3.286; -0.712];
      

 syms x u1 u2 o1 o2;
 
 Pw1=1/3;
 Pw2=2/3;
 
 %Costruisco la Mistura
 Fun_Mixture=@(x,u,o,Pw) Pw(1)*GaussianaMulti(o(1),u(1),x) + Pw(2)*GaussianaMulti(o(2),u(2),x);
 
 n=length(Dati);
 
 %Definisco i settaggi dell' OTTIMIZAZIONE
 options = optimset('Display','final','MaxFunEvals',1500,'MaxIter',1500,'TolFun',1e-10,'TolX',1e-10);
 
 b=max(Dati);
 a=min(Dati);
 %Inizializzo i Valori della MEDIA con 2 numeri estratti a caso nell'
 %intervallo dei campioni
 Uo=a+(b-a)*rand(1,2);
 
 %Inizializzo la Varianza con 2 valori casuali in [0 1]
 Oo=2+randn(1,2);
 
 %Inizializzo le Probabilit� a priori
 Pwo=[.5 .5];
 
 Vo=[Oo Uo Pwo];
 
 %Calcolo il Massimo della Funzione Verosimiglianza per la Mistura
 %Paramtrizzata
 fun=@(X) -Mixture_LogLike(X(1:2).^2,X(3:4),Dati,X(5),X(6));
 
 x0 = [-1,1];     % Make a starting guess at the solution
options = optimset('LargeScale','off','Display','iter','MaxIter',20);
[V, fval] = fmincon(fun,Vo,[],[],[],[],[],[],@CondizioniPriori,options)
 
 
 %[V ,val,exitflag,output]=fminsearch(fun,Vo,options)
 Media1=V(3:4);
 Var1=V(1:2);
 Pw1=V(5:6);
 
 syms x
 F=Fun_Mixture(x,Media1,Var1,Pw1);
 ezplot(F);
 hold on; 
 
 plot(Dati,zeros(1,25),'og','MarkerSize',2,'MarkerFaceColor','g');
 plot(Media1(1),0,'or','MarkerSize',2,'MarkerFaceColor','r')
 plot(Media1(2),0,'or','MarkerSize',2,'MarkerFaceColor','r')
 title('Mistura di Gaussiane con Media Stimata mediante ML- FMINSEARCH')
 
 display('Metodo Iterativo: premere un tasto per continuare..')
 pause;
 %Applico l' algoritmo Iterativo per la Stima del Vettore delle Medie:
 %
 U=Uo.';
 O=(Oo.^2).';
 Pw=Pwo.';
 J=1;
 Condizione=0;
 maxIter=10;
 while (J<maxIter)
     %U(:,J+1)=0;
     Den=0;
     P1=zeros(1,25);
     P2=zeros(1,25);
     for K=1:25 
         
        % Valuto la Mistura nel Punto x(k)  P(x(k) | �)
        p(K)=Mixture(O(:,J),U(:,J),Dati(K),Pw(1,J),Pw(2,J));
        
        
        %Calcolo la Verosimiglianza del Campione Xk alla componente 1 della
        %Mistura: p(w(i) | x(k) , �(i))
        %
        % p(x(k)|w(1) , �(1))
        p1(K)=GaussianaMulti(O(1,J),U(1,J),Dati(K));
        
        %Calcolo la Verosimiglianza del Campione Xk alla componente 2 della
        %Mistura: p(w(i) | x(k) , �(i))
        %
        % p(x(k)|w(1) , �(2))
        p2(K)=GaussianaMulti(O(2,J),U(2,J),Dati(K));
        
        %Calcolo la Verosimiglianza del Campione Xk alla componente k della
        %Mistura ->  P(w(k) | x(k) , �(k))
        P1(K)=p1(K)*Pw(1,J)/p(K);
        P2(K)=p2(K)*Pw(2,J)/p(K);
        
     end
     
     %% M-Step
          
         %Aggiorno la Media 
         U(1,J+1)=(P1*Dati)/sum(P1);
         U(2,J+1)=(P2*Dati)/sum(P2);

         %Aggiorno la Varianza utilizzando la Media Corrente
         O(1,J+1)=(P1*(Dati-U(1,J+1)).^2)/sum(P1);
         O(2,J+1)=(P2*(Dati-U(2,J+1)).^2)/sum(P2);
         
         %Aggiorno i Prior
         Pw(1,J+1)=sum(P1)/25;
         Pw(2,J+1)=sum(P2)/25;

         disp('.')
     
     J=J+1;
 end
 Media2=U(:,end).';
 Var2=O(:,end).';
 Pw2=Pw(:,end).';
 
 figure;
 syms x
 F=Fun_Mixture(x,Media2,Var2,Pw2);
 ezplot(F,[-4 4]);
 hold on;
 plot(Dati,zeros(1,25),'og','MarkerSize',2,'MarkerFaceColor','g');
 plot(Media2(1),0,'or','MarkerSize',2,'MarkerFaceColor','r')
 plot(Media2(2),0,'or','MarkerSize',2,'MarkerFaceColor','r')
 title('Mistura di Gaussiane con Media Stimata mediante ML iterativa')
 axis tight