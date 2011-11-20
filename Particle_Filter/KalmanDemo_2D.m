%% Simulazione di Kalman
%
clc;clear;close all
%% PROCESSO DI STATO
%
% X(k+1)=A*X(k) + w(k)
%
% A:        Matrice di Transizione di Stato
% X(k):     Stato al tempo precedente
% G:        Matrice Evolutiva del Rumore
% w(k):     Rumore di Processo ( Var Gaussiana --> N(0,Q) )
%Matrice di Evoluzione di Stato                      
A = [1 2 0 0;
     0 1 0 0;
     0 0 1 3;
     0 0 0 1];
%Matrice di Covarianza del Rumore di Processo
Q=4*eye(4);

 
%% PROCESSO DI MISURA
%
% Y(k)=H*X(k) + v(k)
%
% H:        Matrice che lega lo STATO alle OSSERVAZIONI
% X(k):     Stato al tempo corrente
% v(k):     Rumore di MISURA ( Var Gaussiana --> N(0,R) )

%Matrice di Relazione tra STATO ed OSSERVAZONE
%nb:
%
% C*X(k)= [X(k,1) X(k,3)] mette a zero le componenti vx e vy
H = [1 0 0 0;
     0 0 1 0]; 
 
%Varianza
target_delta=700; 
%Covarianza del RUMORE di misura
R=[target_delta^2  0;
    0  target_delta^2]; 
                          

%% Genero L'evoluzione del Sistema Dinamico che cerco di ricostruire:
%  1)assegno lo Stato Inziale 
%  2)faccio evolvere il sistema
n=30;
stato=zeros(4,n);   
misura=zeros(2,n);
%STATO Iniziale: [x vx y vy]
stato(:,1)=[1500 500 1500  400].';

for i=2:n
       %Evoluzione dello STATO del sistema Dinamico senza ERRORE 
       stato(:,i)=A*stato(:,i-1);     
       %misure generate dal Modello
       misura(:,i)=H*stato(:,i)+sqrt(R)*(randn(2,1));
end
plot(stato(1,:),stato(3,:),'-');hold on
plot(misura(1,:),misura(2,:),'or');

xlabel('x(k)'),ylabel('y(k)');


s.x =stato(:,1) ;   % Stato Iniziale
s.A = A;            % Matrice di Transizione
% Define a process noise (stdev) of 2 volts as the car operates:
s.Q = Q; % variance, hence stdev^2

% Define the voltimeter to measure the voltage itself:
s.H = H;
% Define a measurement error (stdev) of 2 volts:
s.R = R; % variance, hence stdev^2

% Do not define any system input (control) functions:
s.B = 0;
s.u = 0;

s.P = eye(4);
% Do not specify an initial state:
% Generate random voltages and watch the filter operate.
 
for t=1:n
    
   s(t).z = misura(:,t);   % create a measurement
   s(t+1) = kalmanf(s(t)); % perform a Kalman filter iteration
   statoStimato(:,t)=s(t+1).x;
end
 
% plot measurement data:
 plot(statoStimato(1,:),statoStimato(3,:),'-r');
 legend('Processo da Ricostruire ','Misure','Kalman Estimate',4)

