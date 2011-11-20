% Extended Kalman Filter Demo

%{

Version 1.0, September 2006

This tutorial was written by Jose Manuel Rodriguez, 
Universidad Politecnica de Catalu�a, Spain
Based on the Technical Document by
Juan Andrade-Cetto, "The Kalman Filter", March 2002
%}

%{
x is the real plant behavior, in this case a sinus wave
with the following formulation:

    x(time)=sin(frec(@ time-1)*time-1) + ramdom error signal("sigmav")
    frec(time)=frec(@ time-1)

x_ is the predicted state, this is where Kalman filter will come 
and where we will correct our estimations using an observation

z is the observation of the real plant, in this case corresponding only to 
the position at a given time. Note that this observation is subject to an
error, therefore the resulting equation is: 

        z(time)=x(time)+ramdom error("sigmaw")

Our first prediction will come from the plant ideal behavior, then
using the observations and error covariance we will obtain a better estimate.

xc is the ideal plant behavior... this is used just for comparison

P is the state error covariances at a given time of all the involved variables,
note that we are forming x as a 2 by 1 matrix with the following:
x(1,n) -> position
x(2,n) -> frecuency

Our functions are as following for this example:
Let's say:
f1: x1(time)=sin(x2(time-1)*(time-1))+V     V->ramdom plant error
f2: x2(time)=x2(time-1)
h:  y=x1+w                                  w->ramdom sensor error

F is the Jacobian of the transfer function due to the involved variables,
in this case these are x1 and x2, therefore F will be a 2 by 2 matrix 
(always the matrix is square). The resulting F depends on time and must be
computed for every step that the system takes.

F is as follows:
F -> df1/dx1 = 0    df1/dx2 = cos(x2*time)*time
     df2/dx1 = 0    df2/dx1 = 1

H=dh/dx1 = 1
  dh/dx2 = 0
%}

clear all; close all;


% Initial Conditions
x(:,1)  = [0;0.05];     %Our real plant initial condition       
x_(:,1) = [0;0.04];     %Our estimate initial conidition (they might differ)

xc = x_;                %Set the ideal model as we think it should start
P = [0.01 0;            %set initial error covariance for position & frec, both at sigma 0.1, P=diag([sigma_pos_init^2 sigmav_frec_init^2])
     0     0.01];
sigmav = 0.1;           %the covariance coeficient for the position error, sigma
sigmaw = 0.5;           %the covariance coeficient for the frecuency error, sigma
Q = sigmav*sigmav;      %the error covariance constant to be used, in this case just a escalar unit
R = sigmaw*sigmaw;      %the error covariance constant to be used, in this case just a escalar unit

G = [1;0];              %G is the Jacobian of the plant tranfer functions due to the error.
H = [ 1 0];             %H is the Jacobian of the sensor transfer functions due to the variables involved
W = 1;                  %W is the Jacobian of the sensor transfer functions due to the error.

steps = 1000;   %Amount of steps to simulate

% bucle
for i =2:steps          %start @ time=2 
  % the real plant
  x(:,i) = [sin(x(2,i-1)*(i-1)) + randn*sigmav ; 
            x(2,i-1) ];
        
  z(i) = x(1,i) + randn*sigmaw;
  
  % prediction: da notare l' assenza del RUMORE di STATO e MISURA
  x_(:,i) = [sin(x_(2,i-1)*(i-1)); x_(2,i-1)];
  z_(i) = x_(1,i);

  % compute F
  F = [0 i*cos(x_(2,i)*i);
       0 1];
  
  % Prediction of the plant covariance
  P = F*P*F' + G*Q*G';
  % Innovation Covariance
  S = H*P*H'+R;
  % Kalman's gain
  K = P*H'*inv(S);
  % State check up and update
  x_(:,i) = x_(:,i) + K * (z(i)-z_(i));
  
  % Covariance check up and update
  P = (eye(2)-K*H)*P;
  
  sigmaP(:,i)=sqrt(diag(P)); %sigmap is for storing the current error covariance for ploting pourposes
end

figure(1);clf; hold on;
plot(x(1,:),'-b');                  %plot the real plant behavior
plot(z,'.r');                       %plot the observations over this plant

plot(x_(1,:),'-r');                 %plot the Kalman filter prediction over the plant
 

%These two are the threshold in witch I'm certain that the plant state is at a given time
plot(x_(1,:)+2*sigmaP(1,:),'-g');  
plot(x_(1,:)-2*sigmaP(1,:),'-g');


figure(2);clf;hold on;
plot(x(2,:),'-b');                  %Frecuency estimation
plot(x_(2,:),'-g');                 %Frecuency filtered by Kalman
  