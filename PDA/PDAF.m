%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function  PDAF(target_position,n,T,MC_number)
% Parametri IN :
%           -target_position:  Punto Inziale
%           - n:  Numero di Evoluzioni che desidero ricostruire 
%           - T:  Tempo dello Stato iniziale
%           - MC_number: Numero simulazioni che effettuo
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  PDAF(target_position,n,T,MC_number)
tic
%Probabilità che descrive le Misure
%Probabilità di Misura Corretta
Pd=1;   
%Probabilità che la mis
Pg=0.99;  
% Questa e la Distanza Massima consentita ( g^2)
g_sigma=9.21;                                            
gamma=1*10^(-6);                                         
%Dimensione dello Spazio Misure
c2=pi;
%% PROCESSO DI STATO
%
% X(k+1)=A*X(k) + G*w(k)
%
% A:        Matrice di Transizione di Stato
% X(k):     Stato al tempo precedente
% G:        Matrice Evolutiva del Rumore
% w(k):     Rumore di Processo ( Var Gaussiana --> N(0,Q) )
%Matrice di Evoluzione di Stato                      
A = [1 T 0 0;
     0 1 0 0;
     0 0 1 T;
     0 0 0 1];
%Matrice di Covarianza del Rumore di Processo
Q=[4 0;
   0 4]; 

%Matrice che definisce l'evoluzione del RUMORE nel processo di Evoluzione
%di stato
G=[T^2/2 0;
       T 0;
       0 T^2/2;
       0 T]; 
%% PROCESSO DI MISURA
%
% Y(k)=C*X(k) + v(k)
%
% C:        Matrice che lega lo STATO alle OSSERVAZIONI
% X(k):     Stato al tempo corrente
% v(k):     Rumore di MISURA ( Var Gaussiana --> N(0,R) )

%Matrice di Relazione tra STATO ed OSSERVAZONE
%nb:
%
% C*X(k)= [X(k,1) X(k,3)] mette a zero le componenti vx e vy
C = [1 0 0 0;
     0 0 1 0]; 
 
%Varianza
target_delta=100; 
%Covarianza del RUMORE di misura
R=[target_delta^2  0;
    0  target_delta^2]; 
                                        
% Stima Iniziale della Matrice di Cov Errore di Stima sullo Stato
P=[target_delta^2 target_delta^2 0 0;
   target_delta^2 2*target_delta^2 0 0;
   0 0 target_delta^2 target_delta^2;
   0 0 target_delta^2 2*target_delta^2];
                                                         
%Valori Ricostruiti dal Filtro
x_filter=zeros(4,n);                                     
x_filter1=zeros(4,n,MC_number);                          

%% Genero L'evoluzione del Sistema Dinamico che cerco di ricostruire:
%  1)assegno lo Stato Inziale 
%  2)faccio evolvere il sistema
data_status=zeros(4,n);   

data_status(:,1)=target_position';

for i=2:n
       %Evoluzione dello STATO del sistema Dinamico 
       data_status(:,i)=A*data_status(:,i-1)+G*sqrt(Q)*(randn(2,1));   %ʵ��λ��  
end

figure(2)
plot(data_status(1,:),data_status(3,:),'-o');hold on
xlabel('x(m)'),ylabel('y(m)');
legend('Processo da Ricostruire ',4)
%axis([0 30 1 25])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PROCESSO DI TRACKING:
% Cerco di ricostruire la Traiettoria (x,y) attraverso PDA
%  
% 1) Genero delle Osservazioni Multiple y
% 2) Filtro le Osservazioni y attraverso il GATE gsigma
% 3) Applico il Probabilistic Data Association
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

for M=1:MC_number

    Noise=[];
    data_measurement=zeros(2,n);
    %% Genero le n OSSERVAZIONI secondo il MODELLO di MISURA
    %  solo la Posizione (x,y) essendo:
    %
    % C*X(k) +  v(k) = [X(k,1) X(k,3)] +  v(k)
    for i=1:n        
         data_measurement(1,i)=data_status(1,i)+randn(1)*target_delta;
         data_measurement(2,i)=data_status(3,i)+randn(1)*target_delta;     
    end
    NOISE=[];
    
    %% Genero le Osservazioni MULTIPLE (insieme y) al quale applicare PDA
    for t=1:n
        % Processo di PREDICT and UPDATE del filtraggio Bayesiano
        % Evoluzione di STATO PRIVA DI RUMORE GAUSSIANO
        if t~=1
            %STATO PREDETTO
            x_predic = A*x_filter(:,t-1);                                  %��ǰһʱ�̵��˲�ֵ��Ԥ�⵱ǰ��ֵ 
        else
            %STATO INIZIALE
            x_predic = target_position';                                   %��һ�β�����������ʵλ�õ�Ԥ��ֵ 
        end

        %Matrice Covarianza Errore Predetto
        %P(k|k-1)
        P_predic = A*P*A'+G*Q*G';
        
        %OSSERVAZIONE PREDETTA
        Z_predic = C*x_predic;
        
        %Matrice Covarianza Errore di Misura: S(k)
        S = C*P_predic*C'+ R;
       
        %Matrice di Gain
        K = P_predic*C'*inv(S);                                            %����
        
        %Calcolo VALIDATION GATE
        %
        %Volume VG: siamo nel 2D
        ellipse_Volume=c2*g_sigma*sqrt(det(S));
        disp(det(S))
        inv(S)
        %number_returns=floor(10*ellipse_Volume*gamma+1)
          
        %GENERO UN UPPERBOUND PER IL VALIDATION GATE e campiono attraverso
        %la tecnica di Rigetto
        %Lato del Volume
        side=sqrt((10*ellipse_Volume*gamma+1)/gamma)/2;                   
       
        %% GENERAZIONE Osservazioni MULTIPLE da filtrare attraverso il VALIDATION GATE
        %  Utilizzo MCMC per campionare una Distribuzione CIRCOLARE UNIFORME
        %  centrata attorno alla MISURA PREDETTA o STATO PREDETTO
        %Inizializzo tutti le Variabili presenti nel Modello
        NVar=2;
        Trial=40;
        Pdf=@(Theta)DistrCircolare(Theta,[x_predic(1);x_predic(3) ],40);


        %INIZIALIZZO I CAMPIONI 
        %Campioni che genero nel tempo: ogni colonna rappresenta il campione
        %estratto dal Processo di campionamento
        Theta=zeros(NVar,Trial);
         %Punto Iniziale
        Theta(1,1)=x_predic(1); 
        Theta(2,1)=x_predic(3) ; 
        %Campiono
        Noise=Fun_MetroPolisHastingsSampler_CW(Pdf,Theta,NVar,Trial,[1 20]);
        NOISE=[NOISE Noise];
        number_returns=size(Noise,2);
        
        %Alle Osservazioni RUMOROSE Aggiungo la MISURA CORRETTA
        MisuraCorretta=zeros(1,2);
        MisuraCorretta(1)=data_measurement(1,t);
        MisuraCorretta(2)=data_measurement(2,t);
        %Osservazioni Multiple y1
        y1=[Noise MisuraCorretta'];
        
        %% FILTRAGGIO delle Osservazioni MULTIPLE mediante GATE
        %Osservazioni Multiple y1 filtrate attraverso il GATE g-sigma
        y=[];
        d=[];
        %Numero di Osservazioni Multiple Finali
        m=0;
        for j=1:(number_returns+1)
            %Distanza Osservazione - Osservazione Predetta
            d=y1(:,j)-Z_predic;
            %Calcolo la Distanza di Malanhobis per scartare le OSSERVAZIONI
            % che non cadono NEL GATE (D < g_sigma)
            D=d'*inv(S)*d;
            if D<=g_sigma
                y=[y y1(:,j)];                                              
                m=m+1;                                                      
            end
        end
        DisegnaValidationGate_Predict
        
        %% PDA FILTERING con y OSSERVAZIONI MULTIPLE al tempo t
        %Calcolo Bk
                                          
        if m==0
           %Tutti le osservazioni sono fuori gate:
           % Lo Stato Stimato e uguale allo stato PREDETTO
           x_filter(:,t)= x_predic;
           P=P_predic; 
           %Disegno 
            DisegnaValidationGate_Update
           
        else   
            %STIMA DELLO STATO via PROBABILISTIC DATA ASSOCIATION
            
            %% Calcolo dei coefficienti BeTa(i)
            %
            % BeTa(i) = e(i)/(b + sum(e(i))) i=1..m
            % BeTa(i) = b/(b + sum(e(i))) i=0
            %
            BeTa=zeros(1,m);
            % Calcolo il coefficiente b da cui dipende BeTa(i)
            b=gamma*2*pi*sqrt(det(S))*(1-Pd*Pg)/Pd;
            
            %Calcolo i Coefficienti e(i) da cui dipende BeTa(i)
            e=zeros(1,m);
            residuo=zeros(2,m);
            E=0;
            for i=1:m
                residuo(:,i)=y(:,i)-Z_predic;
                % Distanza tra OSSERVAZIONE e OSSERVAZIONE PREDETTA
                e(i)=(residuo(:,i))'*inv(S)*(residuo(:,i));
                % quanto e VEROSIMILE L-OSSERVAZIONE y all osservazione
                % predetta
                E=E+exp(-e(i)/2);
            end 
            
            %Posso calcolare BeTa
            BeTa0=b/(b+E);
            %Calcolo
            Vk=zeros(2,1);
            %Termine necessarop per calcolare Spread Of innovetion term   
            v1=zeros(2,2);
            for i=1:m
                BeTa(i)=e(i)/(b+E);   
                 
                %Combined Innovation 
                Vk=Vk+BeTa(i)*(residuo(:,i));
                %Calcolo La somma termine Combined Innovation
                v1=v1+BeTa(i)*(residuo(:,i))*(residuo(:,i))';
            end
            %% Update (FILTRAGGIO BAYESIANO RICORSIVO )
            %Stato Stimato= Combinazione STATO PREDETT0 e STATO
            x_filter(:,t)= x_predic + K*Vk;
            
            %Disegno 
            DisegnaValidationGate_Update
            %Calcolo la MATRICE DI COVARIANZA AGGIORNATA con le Misure
            %Correte
            %originale
            Pc=(eye(4)-K*C)*P_predic;  %--> KALMAN FILTER
            %Pc=P_predic-K*S*K';
            
            %Spread Of innovetion term 
            PP=K*(v1-Vk*Vk')*K';
            
            %Matrice di Covarianza a Posteriori sull'errore di stato
            P=BeTa0*P_predic+(1-BeTa0)*Pc+PP;
            
        end

        x_filter1(:,t,M)=x_filter(:,t);

    end
end
toc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Valore Medio
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R=sum(U1,4)/MC_number;  
x_filter=sum(x_filter1,3)/MC_number;                         
figure                                                        
plot(x_filter(1,:),x_filter(3,:),'r-'),hold on
%plot(data_measurement1(1,:),data_measurement1(3,:),'-')
plot(data_measurement(1,:),data_measurement(2,:),'k-')
plot(NOISE(1,:),NOISE(2,:),'.')                               
%axis([0 30 1 25])
xlabel('x(m)'),ylabel('y(m)');
legend('���Ƶ�λ��','������λ��','�Ӳ�λ��',4)


figure                                                       
plot(1:n,x_filter(1,:),'r-'),hold on
%plot(1:n,data_measurement1(1,:),'-')
plot(data_measurement(1,:),'k-')
%axis([0 50 1 30])
xlabel('t(s)'),ylabel('x(m)');
legend('X�������λ��','X�������λ��',4)


figure                                                       
plot(1:n,x_filter(3,:),'r-'),hold on
%plot(1:n,data_measurement1(3,:),'-')
plot(data_measurement(2,:),'k-')
%axis([0 50 1 25])
xlabel('t(s)'),ylabel('y(m)');
legend('Y�������λ��','Y�������λ��',4)

figure
a=zeros(1,n);                                               
for j=1:n
        a(1,j)=sqrt((x_filter(1,j)-data_status(1,j))^2+(x_filter(3,j)-data_status(3,j))^2);
end
plot(1:n,a(1,:),'r:') 
xlabel('t(s)'),ylabel('Ԥ�����(m)');
figure
plot(1:n,x_filter(2,:),'r-'),hold on
plot(1:n,x_filter(4,:),'-');
xlabel('t(s)'),ylabel('�ٶ�(m/s)');
legend('X�����ٶ�','Y�����ٶ�',4)

