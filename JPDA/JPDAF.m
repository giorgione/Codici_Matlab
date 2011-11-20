function  JPDAF(target_position,n,T,MC_number,c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%           -target_position:  Ŀ��ĳ�ʼλ��
%           - n:  Numero misure effettuatte ad ogni passo 
%           - T:  ������
%           -MC_number:�������
%           - c:  Ŀ�����
 
  
Pd=1;                                                                              %������
Pg=0.99;                                                                           %��ȷ��������������ڵø���
g_sigma=9.21;                                                                      %����
lambda=2; 
gamma=lambda*10^(-6); 
%Misure Effettuate sul Target
Target_measurement=zeros(c,2,n); 

%NOISE VECTOR
target_delta=[100 100];                                                            %Ŀ���Ӧ�Ĺ۲��׼��                    
P=zeros(4,4,c);                                                                    %Э�������
P1=[target_delta(1)^2 0 0 0;
                    0 0.01 0 0;
                    0 0 target_delta(1)^2 0;
                     0 0 0 0.01];        %��ʼЭ������� 
P(:,:,1)=P1;
P(:,:,2)=P1;

A = [1 T 0 0;
     0 1 0 0;
     0 0 1 T;
     0 0 0 1];   
 
C = [1 0 0 0;
     0 0 1 0];
 
%Matrice di Covarianza Associata al Processo di Misura 
R=[target_delta(1)^2 0;
    0 target_delta(1)^2]; 
%Matrice di Covarianza Associata al Processo di Evoluzione
Q=[4 0;
   0 4];                                                  
%
G=[T^2/2 0;
       T 0;
       0 T^2/2;
       0 T];
% Matrice c(numero tracce) x 2(numero features) x n (numero misure):
% conterra' lo stato ricostruito
x_filter=zeros(4,c,n); 
% Campioni generati per lo stato corrente
x_filter1=zeros(4,c,n,MC_number);     

% Matrice c(numero tracce) x 2(numero features) x n (numero misure)
data_measurement=zeros(c,2,n);                                                     %�۲�洢����
data_measurement1=zeros(c,4,n);                                                    %ʵ��λ�����x,y����   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GENERO le 2 Tracce che Desidero Ricostruire Facendo Evolvere il Sistema
% dinamico utilizzando il punto iniziale target_position
%
% STATO TRUE
data_measurement1(:,:,1)=target_position;                                          %ʵ��λ�þ����ʼ�� 
%TRACKS 1 - 2
for i=1:c
    %Faccio Evolvero lo STATO DEL SISTEMA
    for ii=2:n                                                                     %ʵ��λ�� 
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
    %Fisso la TRACCIA i-esima ed Eseguo un filtro di Kalman in modo
    %indipendente su di essa
    for i=1:c      
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
        %Matrice di Covarianza sull errore Stimata
        P_predic=A*P(:,:,i)*A'+G*Q*G';                                                  
        
        %Misura stimata
        Z_predic(:,i)=C*x_predic(:,i);                                                  
        
        %
        R=[target_delta(i)^2 0; 
            0 target_delta(i)^2];
        %Stato Stimato
        S(:,:,i)=C*P_predic*C'+R;                                                       
        
        ellipse_Volume(i)=pi*g_sigma*sqrt(det(S(:,:,i)));  
        
        number_returns=floor(ellipse_Volume(i)*gamma+1);                               
        side=sqrt((ellipse_Volume(i)*gamma+1)/gamma)/2;                                 
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
        b=zeros(1,2);
        b(1)=data_measurement(1,1,t);   %x
        b(2)=data_measurement(1,2,t);   %y
        y1=[NOISE b'];                                                                 %�����յ������еĻز�����y1��,�����Ӳ��͹۲�
        b(1)=data_measurement(2,1,t);
        b(2)=data_measurement(2,2,t);
        y1=[y1 b'];                                                                    %����һ���Ӳ��ز�ʱ
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  3. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        m1=0;                                                                          %��¼��Ч�۲����
        [n1,n2]=size(y1);
        Q1=zeros(100,3);                                                               %Q1Ϊ�۾���
        for j=1:n2 
            flag=0;
            for i=1:c
                d=y1(:,j)-Z_predic(:,i);                                              %�����лز�λ��-Ԥ��λ�ã�����в�����
                D=d'*inv(S(:,:,i))*d;                                                 %��в������ķ���       
                if D<=g_sigma                                                         %������ڡ���                                                   
                   flag=1;
                   Q1(m1+1,1)=1;                                                      %�ӵ�һ�п�ʼ��Q1��ֵ                                                                                            
                   Q1(m1+1,i+1)=1;
                end
            end
                if flag==1   
                   y=[y y1(:,j)];                                                      %������������е����лز�����y��
                   m1=m1+1;                                                            %��¼��Ч�۲����
                end
        end
        Q2=Q1(1:m1,1:3);                                                               %Q2�洢��Ч�ز�
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  4.����������A_matrix,����num��ʾ���������¼�����  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A_matrix=zeros(m1,3,10000);
        A_matrix(:,1,1:10000)=1;
        if m1~=0                                                                       %m1=0��ʾ����Ŀ�궼û�й۲�
           num=1;

           for i=1:m1
            if Q2(i,2)==1                                                              %�ж�������aĿ��Ļز�
                A_matrix(i,2,num)=1;A_matrix(i,1,num)=0;                               %num=1��¼������aĿ��Ļز����ĸ�
                num=num+1;
                for j=1:m1                                                              
                    if (i~=j)&(Q2(j,3)==1)
                        A_matrix(i,2,num)=1;A_matrix(i,1,num)=0;
                        A_matrix(j,3,num)=1;A_matrix(j,1,num)=0;
                        num=num+1;
                    end
                end
            end
        end                                   

        for i=1:m1
            if Q2(i,3)==1
                A_matrix(i,3,num)=1;A_matrix(i,1,num)=0;
                num=num+1;
            end
        end

        else
            flag=1;
        end
        A_matrix=A_matrix(:,:,1:num);                                                  %��ٷ���ֵĽ�����A_matrix��
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  5.����������Pr,����False_num��ʾ������,mea_indicator��ʾ�۲�ָʾ��,target_indicator��ʾĿ��ָʾ��  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pr=zeros(1,num);
    for i=1:num
        False_num=m1;
        N=1;
        for j=1:m1
            mea_indicator=sum(A_matrix(j,2:3,i));                                      %�ο�������ʽ4-48
            if mea_indicator==1                                                        %
                False_num=False_num-1;
                if A_matrix(j,2,i)==1                                                  %���۲���Ŀ��1����
                    b=(y(:,j)-Z_predic(:,1))'*inv(S(:,:,1))*(y(:,j)-Z_predic(:,1));    %��в��
                    N=N/sqrt(det(2*pi*S(:,:,1)))*exp(-1/2*b);                          %������̬�ֲ�����                         
                else                                                                   %���۲���Ŀ��2����
                    b=(y(:,j)-Z_predic(:,2))'*inv(S(:,:,2))*(y(:,j)-Z_predic(:,2));
                    N=N/sqrt(det(2*pi*S(:,:,2)))*exp(-1/2*b);                          %������̬�ֲ�����                         
                end                                                                        
            end
        end
        if Pd==1
            a=1;
        else
            a=1;
            for j=1:c
                target_indicator=sum(A_matrix(:,j+1,i));                               %�ο�������ʽ4-49
                a=a*Pd^target_indicator*(1-Pd)^(1-target_indicator);                   %���������
            end
        end                                                                            
        V=ellipse_Volume(1)+ellipse_Volume(2);                                         %��ʾ�����������

        a1=1;
        for j=1:False_num
            a1=a1*j;
        end
        Pr(i)=N*a*a1/(V^False_num);
    end
    Pr=Pr/sum(Pr);
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  6.�����������U  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    U=zeros(m1+1,c);
    for i=1:c
        for j=1:m1
            for k=1:num
                U(j,i)=U(j,i)+Pr(k)*A_matrix(j,i+1,k);
            end
        end
    end
    U(m1+1,:)=1-sum(U(1:m1,1:c));                                                      %��������Ŀ��T�����Ĺ������ʴ���U��m1+1,:),��һ��
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  7.Kalman�˲���ʼ  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:c                                                                          %����Э�������
        P_predic = A*P(:,:,i)*A'+G*Q*G';
        K (:,:,i)= P_predic*C'*inv(S(:,:,i));
        P(:,:,i)= P_predic-(1-U(m1+1,i))*K(:,:,i)*S(:,:,i)*K(:,:,i)';
    end
    for i=1:c
        a=0;         
        b=0;
        x_filter2=0;                                                                   %������õ��м����
        for j=1:m1
            x_filter2=x_filter2+U(j,i)*(x_predic(:,i)+ K (:,:,i)*(y(:,j)- Z_predic(:,i)));
        end
        x_filter2=U(j+1,i)*x_predic(:,i)+x_filter2;
        x_filter(:,i,t)=x_filter2;
        for j=1:m1+1
            if j==m1+1
                a=x_predic(:,i);
            else
               a=x_predic(:,i)+ K (:,:,i)*(y(:,j)- Z_predic(:,i));
            end
            b=b+U(j,i)*(a*a'-x_filter2*x_filter2');
        end
        P(:,:,i)=P(:,:,i)+b; 
        x_filter1(:,i,t,M)=x_filter(:,i,t);   
    end
    end
    end

    %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%
    %%%%%  ��ͼ  %%%%%
    %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%
    x_filter=sum(x_filter1,4)/MC_number;                                               %�˲�ֵ��ƽ��
    %%%%%%%%%%%%%%%%%%%%
    %%%  1.�˲����  %%%
    %%%%%%%%%%%%%%%%%%%%
    figure;
    %Ŀ��a,b�Ĺ۲�λ��
    for i=1:c
        a=zeros(1,n);
        b=zeros(1,n);
        for j=1:n
            a(j)=data_measurement(i,1,j);
            b(j)=data_measurement(i,2,j);
        end
        if i==1
           plot(a(:),b(:),'b:')
        else 
           plot(a(:),b(:),'b:')
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
legend('Ŀ��a�Ĺ۲�λ��','Ŀ��b�Ĺ۲�λ��','Ŀ��a���Ӳ�','Ŀ��b���Ӳ�','Ŀ��a�Ĺ���λ��','Ŀ��b�Ĺ���λ��',4);grid;
%%%%%%%%%%%%%%%%%%%%
%%%  2.�ٶ����  %%%
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
temp=c1(1,:);
a_extra=zeros(2,n);
b_extra=zeros(1,n);
c_extra=zeros(1,n);
a_extra(1,:)=temp;
a_extra(2,:)=1:1:n;
b_extra=a_extra(1,:);
[c_extra,pos]=sort(b_extra);                                                       %posΪ�������±�,cΪ��һ�е�������;
a_extra(2,:)=a_extra(2,pos);                                                       %�ڶ��а��յ�һ��������±��Ӧ
a_extra(1,:)=c_extra;                                                              %��һ�н�����¸���a �ĵ�һ��;
str1=num2str(a_extra(2,n));
str2=num2str(a_extra(1,n));
str=strcat('\itN=',str1,'\itError=',str2,'(m)');
text(a_extra(2,n),0.8*a_extra(1,n),str);
hold on;
plot([a_extra(2,n) a_extra(2,n)],[0 a_extra(1,n)],'r');
hold on;
plot(1:n,c1(1,:),'r:') 
hold on;
a=0;
for j=1:n
    for i=1:MC_number                                                              %��С�����
        a=(x_filter1(1,2,j,i)-data_measurement1(2,1,j))^2+(x_filter1(3,2,j,i)-data_measurement1(2,3,j))^2;
        c1(2,j)=c1(2,j)+a;
    end
        c1(2,j)=sqrt(c1(2,j)/MC_number);
end
temp=c1(2,:);
a_extra=zeros(2,n);
b_extra=zeros(1,n);
c_extra=zeros(1,n);
a_extra(1,:)=temp;
a_extra(2,:)=1:1:n;
b_extra=a_extra(1,:);
[c_extra,pos]=sort(b_extra);                                                       %posΪ�������±�,cΪ��һ�е�������;
a_extra(2,:)=a_extra(2,pos);                                                       %�ڶ��а��յ�һ��������±��Ӧ
a_extra(1,:)=c_extra;                                                              %��һ�н�����¸���a �ĵ�һ��;
str1=num2str(a_extra(2,n));
str2=num2str(a_extra(1,n));
str=strcat('\itN=',str1,'\itError=',str2,'(m)');
text(a_extra(2,n),0.8*a_extra(1,n),str);
hold on;
plot([a_extra(2,n) a_extra(2,n)],[0 a_extra(1,n)],'b');
hold on;
plot(1:n,c1(2,:),'b:') 
xlabel('times'),ylabel('����ֵ�����ֵ���/m');
legend('Ŀ��a��������ֵ','Ŀ��a�����','Ŀ��b��������ֵ','Ŀ��b�����');grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Revised on 26th June 2008 by wangzexun  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
extra11=zeros(1,n);
extra12=zeros(1,n);
extra13=zeros(1,n);
for j=1:n
    extra11(1,j)=sqrt(x_filter(1,1,j)-data_measurement1(1,1,j))^2+(x_filter(3,1,j)-data_measurement1(1,3,j))^2;
    extra12(1,j)=sqrt((data_measurement(1,1,j)-data_measurement1(1,1,j))^2+(data_measurement(1,2,j)-data_measurement1(1,3,j))^2);
    extra13(1,j)=extra12(1,j)/extra11(1,j);
end
plot(1:n,extra13(1,:),'k:'); 
xlabel('times'),ylabel('RMSE of a');
grid;

figure;
extra21=zeros(1,n);
extra22=zeros(1,n);
extra23=zeros(1,n);
for j=1:n
    extra21(1,j)=sqrt(x_filter(1,2,j)-data_measurement1(2,1,j))^2+(x_filter(3,2,j)-data_measurement1(2,3,j))^2;
    extra22(1,j)=sqrt((data_measurement(2,1,j)-data_measurement1(2,1,j))^2+(data_measurement(2,2,j)-data_measurement1(2,3,j))^2);
    extra23(1,j)=extra22(1,j)/extra21(1,j);
end
plot(1:n,extra23(1,:),'k:'); 
xlabel('times'),ylabel('RMSE of b');
grid;