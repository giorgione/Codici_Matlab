function  JPDAF(target_position,n,T,MC_number,c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function  PDAF(target_position,n,T,MC_number)
% Parametri IN :
%           -target_position:  Punto Inziale
%           - n:  Numero di Evoluzioni che desidero ricostruire 
%           - T:  Tempo dello Stato iniziale
%           - MC_number: Numero simulazioni che effettuo
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
Pd=1;                                                                               
Pg=0.99;                                                                            
g_sigma=9.21;                                                                     
lambda=2; 
gamma=lambda*10^(-6); 

%Misure Effettuate sul Target
Target_measurement=zeros(c,2,n); 

%NOISE VECTOR
target_delta=[100 100];                                                                                
P=zeros(4,4,c);                                                                    
P1=[target_delta(1)^2 0 0 0;
                    0 0.01 0 0;
                    0 0 target_delta(1)^2 0;
                     0 0 0 0.01];        
P(:,:,1)=P1;
P(:,:,2)=P1;

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
 
%Matrice di Covarianza Associata al Processo di Evoluzione
Q=[4 0;
   0 4];                                                  

% Matrice assciata che determina l'evoluzione del Rumore nel Tempo
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
%Matrice di Covarianza Associata al Processo di Misura 
R=[target_delta(1)^2 0;
    0 target_delta(1)^2]; 

% Matrice c(numero tracce) x 2(numero features) x n (numero misure):
% conterra' lo stato ricostruito
x_filter=zeros(4,c,n); 
% Campioni generati per lo stato corrente
x_filter1=zeros(4,c,n,MC_number);     

% Matrice c(numero tracce) x 2(numero features) x n (numero misure)
data_measurement=zeros(c,2,n);                                                     
data_measurement1=zeros(c,4,n);                                                    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GENERO le 2 Tracce che Desidero Ricostruire Facendo Evolvere il Sistema
% dinamico utilizzando il punto iniziale target_position
%
% STATO TRUE
data_measurement1(:,:,1)=target_position;                                           
%TRACKS 1 - 2
for i=1:c
    %Faccio Evolvero lo STATO DEL SISTEMA
    for ii=2:n                                                                      
        data_measurement1(i,:,ii)=(A*data_measurement1(i,:,ii-1)')'+(G*sqrt(Q)*(randn(2,1)))'; 
        
    end
end
%Disegno la TRACK 1
a=zeros(1,n);
b=zeros(1,n);
for i=1:n
    a(i)=data_measurement1(1,1,i);
    b(i)=data_measurement1(1,3,i);
end
plot(a,b,'b-');
hold on;

%Disegno la TRACK 2
a=zeros(1,n);
b=zeros(1,n);
for i=1:n
    a(i)=data_measurement1(2,1,i);
    b(i)=data_measurement1(2,3,i);
end
plot(a,b,'r-');
xlabel('x(m)'),ylabel('y(m)');
legend('Track 1','Track2',1);
grid;

%%%%%%%%%%%%%%%%%%%%%%
% PROBLEMA DI MULTI TRACKING 
for M=1:MC_number
    %%%%%%%%%%%%%%%%%%%%    
    %%%  1.  PROCESSO DI MISURA
    %
    %
    %%%%%%%%%%%%%%%%%%%%
    Noise=[];
    %Genero 50 misurazioni a partire dallo  STATO TRUE attraverso il processo 
    %di misurazione 
    for i=1:n
        %Fisso la TRACCIA j-esima e genero le osserevazioni
        for j=1:c                                                                     
            data_measurement(j,1,i)=data_measurement1(j,1,i)+rand(1)*target_delta(j);
            data_measurement(j,2,i)=data_measurement1(j,3,i)+rand(1)*target_delta(j); 
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  2. Attraverso le misure posso ricostruire le TRACKS ORIGINARIE
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S=zeros(2,2,c);
    Z_predic=zeros(2,2);                                                               %�洢����Ŀ��Ĺ۲�Ԥ��ֵ,��ֻ����x,y���
    x_predic=zeros(4,2);                                                               %�洢����Ŀ���״̬Ԥ��ֵ,������x,y����x,y�����ٶ�
    ellipse_Volume=zeros(1,2);
    NOISE_sum_a=[];                                                                    %�洢Ŀ��1���Ӳ�
    NOISE_sum_b=[];  

for t=1:n
    y1=[];
    y=[];
    Noise=[];
    NOISE=[];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Fisso la TRACK i-esima ed Eseguo  il PROCESSO DI TRACKING:
    % Cerco di ricostruire la Traiettoria (x,y) attraverso PDA
    %  
    % 1) Genero delle Osservazioni Multiple y
    % 2) Filtro le Osservazioni y attraverso il GATE gsigma
    % 3) Applico il Probabilistic Data Association
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    for i=1:c      
        %% Genero le Osservazioni MULTIPLE (insieme y) al quale applicare PDA
        % Processo di PREDICT and UPDATE del filtraggio Bayesiano
        % Evoluzione di STATO PRIVA DI RUMORE GAUSSIANO
        if t~=1
            % EVOLUZIONE DELLO STATO ATTRAVERSO IL SISTEMA
            %Predizione dello Stato
            x_predic(:,i) = A*x_filter(:,i,t-1);                                        
        else
            %LA POSIZIONE INIZIALE E` LO STATO REALE al passo 1
            %Predizione dello Stato
            x_predic(:,i)=target_position(i,:)';                                        
        end
        %% FILTRO DI KALMAN
        %   Matrice Covarianza Errore Predetto
        %   P(k|k-1)
        P_predic=A*P(:,:,i)*A'+G*Q*G';                                                  
        
        %Misura Predetta
        Z_predic(:,i)=C*x_predic(:,i);                                                  
        
        % Matrice di Covarianza Errore che varia nel tempo
        R=[target_delta(i)^2 0; 
            0 target_delta(i)^2];
        
        %Matrice Covarianza Errore di Misura: S(k)
        S(:,:,i)=C*P_predic*C'+R;                                                       
        
        %% GENERAZIONE Osservazioni MULTIPLE da filtrare attraverso il VALIDATION GATE
         
        %Calcolo VALIDATION GATE
        %
        %Volume VG: siamo nel 2D
        ellipse_Volume(i)=pi*g_sigma*sqrt(det(S(:,:,i)));  
        number_returns=floor(ellipse_Volume(i)*gamma+1);                               
        side=sqrt((ellipse_Volume(i)*gamma+1)/gamma)/2;   
        
        %Aggiungo del RUMORE alle predizioni:
        Noise_x=x_predic(1,i)+side-2*rand(1,number_returns)*side;                      
        Noise_y=x_predic(3,i)+side-2*rand(1,number_returns)*side;    
        Noise=[Noise_x;Noise_y];
        NOISE=[NOISE Noise];
        
        if i==1
            NOISE_sum_a=[NOISE_sum_a Noise];
        else
            NOISE_sum_b=[NOISE_sum_b Noise];
        end
    end
   
        %Processo le 2 TRACK Congiuntamente
        b=zeros(1,2);
        b(1)=data_measurement(1,1,t);   %x
        b(2)=data_measurement(1,2,t);   %y
        %Aggiungo la Misura Corretta per la Track 1
        y1=[NOISE b'];    
        %Aggioungo la Misura Corretta per la Track 2
        b(1)=data_measurement(2,1,t);
        b(2)=data_measurement(2,2,t);
        y1=[y1 b']; 
        %����һ���Ӳ��ز�ʱ
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  Estrazione delle MISURE nel Validation GATE: MISURE
        %%%  AMMISSIBILI per le 2 track presenti nell'esempio
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        m1=0;                                                                          %��¼��Ч�۲����
        [n1,n2]=size(y1);
        Q1=zeros(100,3); 
        %Numero totale di Osservazioni
        for j=1:n2 
            flag=0;
            %Considero la TRACK i-esima e la Misura j-esima e
            %verifico che  Zj sia nel Validation gate di Xi -> VG( Xi )
            for i=1:c
                d=y1(:,j)-Z_predic(:,i);                                              %�����лز�λ��-Ԥ��λ�ã�����в�����
                D=d'*inv(S(:,:,i))*d;                                                 %��в������ķ���       
                
                if D<=g_sigma                                                         %������ڡ���                                                   
                   flag=1;
                   %l'osservazione m1 può essere generata dal CLUTTER
                   Q1(m1+1,1)=1;   
                   %l'osservazione m1 può essere generata dal TRACK i
                   Q1(m1+1,i+1)=1;
                end
                
            end
            
            %Se Zj € VG(X1) or VG(X2)
            if flag==1   
               y=[y y1(:,j)]; 
               %m1 è l'indice delle Osservazioni che sono nel Validation
               %Gate e vengono inserite iterativamente in y
               m1=m1+1;                                                            
            end
        end
        
        %Q2 conterrà solo le MISURE AMMISSIBILI: 
        % MATRICE m1 x 3 --> MISURE sulle righe - TARGET sulle colonne
        Q2=Q1(1:m1,1:3);                                                               
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Generazione delle POSSIBILI MATRICI DI ASSOCIAZIONE a partire dalle
    %   MISURE AMMISIBILI
    %
    %       num: Numero totale di MATRICI       
    %
    %       A_matrix: Matrice (m1 x 3 x 100) di tutte le possibili 
    %                 Associazioni dove:
    %                 m1 -> numero di MISURE AMMISSIBILI
    %                 3  -> numero tot. di track (1:clutter) (2-3: tracks)
    %                 10000 -> numero massimo ASSOCIAZIONI AMMISSIBILI
    %                         generate dalla MATRICE Q2 secondo la
    %                         procedura descritta nell'articolo
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A_matrix=zeros(m1,3,10000);
        %Tutte le misure possono essere Generate dalla Track Clutter
        A_matrix(:,1,1:10000)=1;
        
        %Se Esistono MISURE AMMISSIBILI
        if m1~=0 
           %contatore delle Possibili Associazioni
           num=1;
           % PROCEDURA DI GENERAZIONE MATRICI ASSOCIAZIONE:
           %
           % 1) LETTURA PER RIGHE ed ESTRAZIONE DEL PRIM0 ELEMENTO 1
           for i=1:m1
               % Misura i-esima e ammissibile per il TARGET 2
               if Q2(i,2)==1                                                              
                    A_matrix(i,2,num)=1;
                    A_matrix(i,1,num)=0;                                
                    num=num+1;
                    %Lettura
                    for j=1:m1                                                              
                        if (i~=j) & (Q2(j,3)==1)
                            A_matrix(i,2,num)=1;
                            A_matrix(i,1,num)=0;
                            A_matrix(j,3,num)=1;
                            A_matrix(j,1,num)=0;
                            num=num+1;
                        end
                    end
                end
            end                                   

            for i=1:m1
                if Q2(i,3)==1
                    A_matrix(i,3,num)=1;
                    A_matrix(i,1,num)=0;
                    num=num+1;
                end
            end

        else
            %Nessuna MISURA E' AMMISSIBILE
            flag=1;
        end
        %Ho Generato Tutte le POSSIBILI ASSOCIAZIONI
        A_matrix=A_matrix(:,:,1:num);                                                  %��ٷ���ֵĽ�����A_matrix��
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Per ogni EVENTO DI ASSOCIAZIONE calcolo:
    %       
    %            k         k-1 
    %       P(  Y  | X ,  Y    ): LIKELYHOOD delle MISURE OSSERVATE
    %
    %      
    %                     k-1 
    %       P( X | m  ,  Y    ):  PRIOR delle MISURE OSSERVATE
    %
    %    e stimo A posteriori: 
    %
    %                    k  
    %       P(  X  |   Y    ): STIMA A POSTERIORI DELLO STATO
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pr=zeros(1,num);
    %Considero tutte le Possibili MATRICI DI ASSOCIAZIONE
    for i=1:num
        
        %Tutte le Misure possono Essere falsi allarmi 
        False_num=m1;
        N=1;
        %% Fissato EVENTO DI ASSOCIAZIONE e la relativa MATRICE DI ASSOCIAZIONE   
        %  Per ogni Misura j vado a calcolare:
        %
        %       - MEASURAMENT ASSOCIATION INDICATOR - tau(j,X)
        %
        %       - TARGET DETECTION INDICATOR - lambda(t,X)
        %
        % 
        for j=1:m1
            %Calcolo il MEASURAMENT ASSOCIATION INDICATOR per l'evento di
            %Associazione Xi (fisso j e sommo su tutti le Track(2-3))
            mea_indicator=sum(A_matrix(j,2:3,i));                                      %�ο�������ʽ4-48
            
            if mea_indicator==1   
                %Considero a quale TRACK LA MISURA E' STATA ASSOCIATA
                
                %Aggiorno FM(X):  il numero Totale di Misure Errate in X
                False_num=False_num-1;
                
                %Individuo a quale delle 2 TRACKS la misura j e ASSOCIATA
                %rispetto all' ASSOCIAZIONE PLAUSIBILE X e calcolo:
                %
                %                  k-1 
                % P(  yj | Xjt ,  Y    )
                if A_matrix(j,2,i)==1   % Xjt                                               %���۲���Ŀ��1����
                    b=(y(:,j)-Z_predic(:,1))'*inv(S(:,:,1))*(y(:,j)-Z_predic(:,1));    %��в��
                    N=N/sqrt(det(2*pi*S(:,:,1)))*exp(-1/2*b);                          %������̬�ֲ�����                         
                else                                                                   %���۲���Ŀ��2����
                    b=(y(:,j)-Z_predic(:,2))'*inv(S(:,:,2))*(y(:,j)-Z_predic(:,2));
                    N=N/sqrt(det(2*pi*S(:,:,2)))*exp(-1/2*b);                          %������̬�ֲ�����                         
                end                                                                        
            end
        end 
        %% N conterra la Joint Probability
        %
        %      k         k-1 
        % P(  Y  | X ,  Y    ) = 
        % 
        %     ____
        %      ||                      k-1 
        %      ||     P(  yj | Xjt ,  Y    )
        %Calcolo il Vulume totale del Validation GAte
        V=ellipse_Volume(1)+ellipse_Volume(2); 
        LikelyHood=N/(V^False_num);
        
        %% Calcolo la PRIOR
         %                     k-1 
         %       P( X | m  ,  Y    ):  PRIOR delle MISURE OSSERVATE
         %
        if Pd==1
            Prior=1;
        else
            Prior=1;
           
            
            %Calcolo il TARGET ASSOCIATION INDICATOR per l'evento di
            %Associazione Xi (fisso j e sommo su tutti le Track(2-3))
            for j=1:c
                %TARGET INDICATOR PARZIALE: sommo solo 
                target_indicator=sum(A_matrix(:,j+1,i));                               %�ο�������ʽ4-49
                
                Prior=Prior*(Pd^target_indicator)*(1-Pd)^(1-target_indicator);                   %���������
            end
        end 
                                               %��ʾ�����������

        %Calcolo il Numero di Eventi in X per i quali lo stesso set di  
        % TARGETS e' rilevato:
        %
        %   m           m!
        % P         =  ----               
        %  m-FM(X)     FM(m)!
        %
        % Calcolo   FM(m)!
        a1=1;
        for j=1:False_num
            a1=a1*j;
        end
        
        %                    k  
        %       P(  X  |   Y    ): STIMA A POSTERIORI DELLO STATO
        Pr(i)=a1*LikelyHood*Prior;
    end
    
    Pr=Pr/sum(Pr);
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  6. Calcolo dei Coefficienti Beta(j,t)
    %%%%%%%%%%%%%%%%%%%%%%%%%
    BeTa=zeros(m1+1,c);
    for i=1:c
        for j=1:m1
            for k=1:num
                BeTa(j,i)=BeTa(j,i)+Pr(k)*A_matrix(j,i+1,k);
            end
        end
    end
    %Calcolo di Beta(0,t) --> ultima posizione di Beta
    BeTa(m1+1,:)=1-sum(BeTa(1:m1,1:c));                                                      
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  7.Kalman PREDICT   per i 2 FILTRI
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:c         
        %Matrice Covarianza Errore Predetto
        %P(k|k-1)
        P_predic = A*P(:,:,i)*A'+G*Q*G';
        %Matrice di GAIN
        K(:,:,i)= P_predic*C'*inv(S(:,:,i));
        
        %Update della MAtrice di Covarianza: manca il Termine Pk che dipende
        % dalla combined innovation
        P(:,:,i)= P_predic-(1-BeTa(m1+1,i))*K(:,:,i)*S(:,:,i)*K(:,:,i)';
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  7.Kalman   UPDATE per i 2 FILTRI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:c
        a=0;         
        Pk=0;
        
        %% Combined Innovation
        % La formula e' in versione compatta rispetto all'articolo in quanto
        % la matrice di Gain K e' moltiplicata per
        x_filter2=0;    
        
        for j=1:m1
            %Stima a Priori X(k|k)
            a=x_predic(:,i)+ K(:,:,i)*(y(:,j)- Z_predic(:,i));
            
            x_filter2=x_filter2+BeTa(j,i)*(a);
        end
        %Aggiungo il Fattore BeTa0
        x_filter2=BeTa(j+1,i)*x_predic(:,i)+x_filter2;
        x_filter(:,i,t)=x_filter2;
        
        a=0;
        for j=1:m1+1
            if j==m1+1
                %Caso di Associazione a Clutter
                a=x_predic(:,i);
            else
                %Caso di Associazione a TRACK
                
               a=x_predic(:,i)+ K (:,:,i)*(y(:,j)- Z_predic(:,i));
            end
            
            Pk=Pk+BeTa(j,i)*(a*a'-x_filter2*x_filter2');
        end
        
        %Update della MAtrice di Covarianza: manca Pk
        P(:,:,i)=P(:,:,i)+b; 
        
        x_filter1(:,i,t,M)=x_filter(:,i,t);   
    end
end % tempo

end %Simulazioni

    %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%
    %%%%%  Valore Medio sulle simulazioni
    %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%
    x_filter=sum(x_filter1,4)/MC_number;                                               %�˲�ֵ��ƽ��
    %%%%%%%%%%%%%%%%%%%%
    %%%  1.Visualizzazione
    %%%%%%%%%%%%%%%%%%%%
    figure;
    % 
    for i=1:c
        a=zeros(1,n);
        b=zeros(1,n);
        for j=1:n
            a(j)=data_measurement(i,1,j);
            b(j)=data_measurement(i,2,j);
        end
        if i==1
           plot(a(:),b(:),'bo')
        else 
           plot(a(:),b(:),'bo')
        end
        hold on;
end
%Ŀ��a,b���Ӳ�λ��
for i=1:c
    if i==1
       plot(NOISE_sum_a(1,:),NOISE_sum_a(2,:),'c.');
    else
       plot(NOISE_sum_b(1,:),NOISE_sum_b(2,:),'c.');
   end
end
hold on;
%Ŀ��a,b�Ĺ���λ��
for i=1:c
    a=zeros(1,n);
    b=zeros(1,n);
    for j=1:n
        a(j)=x_filter(1,i,j);
        b(j)=x_filter(3,i,j);
    end
    if i==1
        plot(a(:),b(:),'r-');
    else 
        plot(a(:),b(:),'r-');
    end
hold on;
end
xlabel('x/m'),ylabel('y/m');
legend('Track1','Track2','Misure','Misure','Track1 Ricostruita','Track1 Ricostruita',4);grid;
%%%%%%%%%%%%%%%%%%%%
%%%  Visualizzo l'Errore Medio Tra Stima dello Stato Prodotto dal Sistema e
%%%  STATO del Sistema ESATTO (che cerco di ricostruire)
%%%%%%%%%%%%%%%%%%%%
figure;
a=0;
c1=zeros(c,n);
for j=1:n
    for i=1:MC_number                                                              %��С�����
        a=(x_filter1(1,1,j,i)-data_measurement1(1,1,j))^2+(x_filter1(3,1,j,i)-data_measurement1(1,3,j))^2;
        c1(1,j)=c1(1,j)+a;
    end
        c1(1,j)=sqrt(c1(1,j)/MC_number);
end
 
plot(1:n,c1(1,:),'r-') 

hold on;
a=0;
for j=1:n
    for i=1:MC_number                                                              %��С�����
        a=(x_filter1(1,2,j,i)-data_measurement1(2,1,j))^2+(x_filter1(3,2,j,i)-data_measurement1(2,3,j))^2;
        c1(2,j)=c1(2,j)+a;
    end
        c1(2,j)=sqrt(c1(2,j)/MC_number);
end
 
plot(1:n,c1(2,:),'b-') 
xlabel('X_{est}(t)'),ylabel('ErroreMedio(X_{est}(t)-X(t))');
legend('Errore Ricostruzione Track1','Errore Ricostruzione Track2'  );
grid;

 