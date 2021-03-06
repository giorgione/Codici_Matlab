%% Chapter 2. (Example 1)
% Metrppolis per Distribuzioni MultiVariate:
% Campionamento di una Bivariate Exponential Distribution caratterizzata da:
%
% P(theta1,theta2| lambda1,lambda2,lambda) =  
%
%      -(lambda1+lambda)*theta1-(lambda2+lambda)*theta2-lambda*max(lambda1,lambda2) 
%    e
%
% con:
%       theta1,theta2 in [0,8]    
%       lambda1 = 0.5;      lambda2 = 0.1;      lambda = 0.01;
%
% Utilizzo come PROPOSAL:  
%                  t       t-1                 
%           Q(theta | theta   ) =  Distribuzione Uniforme in [0 8] <->
%                                        INDIPENDECE SAMPLER
%
% 
% Approccio al Blockwise updating. 
clc;clear;close all
%% Initialize the Metropolis sampler
T= 5000; % Set the maximum number of iterations
thetamin = [ 0 0 ]; % define minimum for theta1 and theta2
thetamax = [ 8 8 ]; % define maximum for theta1 and theta2
seed=1; rand( 'state' , seed ); randn('state',seed ); % set the random seed
theta = zeros( 2 , T );
% Init storage space for our samples: ASSUMO CHE I CAMPIONI SIANO
% INDIPENDENTI
theta(1,1) = unifrnd( thetamin(1) , thetamax(1) ); % Start value for theta1
theta(2,1) = unifrnd( thetamin(2) , thetamax(2) ); % Start value for theta2

%% Start BlockWise sampling 
t = 1;
while t < T % Iterate until we have T samples
    t = t + 1;
    % Propose a new value for theta --> Variabile MultiDimensionale
    thetastar = unifrnd( thetamin , thetamax );
    
    % Calculate the ACCEPTANCE RATIO (essendo la proposal UNIFORME va via 
    % nel computo del Ratio)
    pratio = bivexp( thetastar(1) , thetastar(2) ) / ...
    bivexp( theta(1,t-1), theta(2,t-1) );

    %Probabilità di Accettazione
    alpha = min( [ 1 pratio ] ); 
    
    % Draw a UNIFORM DEVIATE from [ 0 1 ]
    u = rand; 
    
    if u < alpha % Do we accept this proposal?
        theta(:,t) = thetastar; % proposal becomes new value for theta
    else
        theta(:,t) = theta(:,t-1); % copy old value of theta
    end
end

%% Display histogram of our samples
figure( 1 ); clf;
subplot( 1,2,1 );
nbins = 10;
thetabins1 = linspace( thetamin(1) , thetamax(1) , nbins );
thetabins2 = linspace( thetamin(2) , thetamax(2) , nbins );
hist3( theta' , 'Edges' , {thetabins1 thetabins2} );
xlabel( '\theta 1' ); ylabel('\theta 2' ); zlabel( 'counts' );
az = 61; el = 30;
view(az, el);
title('Sampling ')

%% Plot the theoretical density
subplot(1,2,2);
nbins = 20;
thetabins1 = linspace( thetamin(1) , thetamax(1) , nbins );
thetabins2 = linspace( thetamin(2) , thetamax(2) , nbins );
[ theta1grid , theta2grid ] = meshgrid( thetabins1 , thetabins2 );
ygrid = bivexp( theta1grid , theta2grid );
mesh( theta1grid , theta2grid , ygrid );
xlabel( '\theta 1' ); ylabel('\theta 2' );
zlabel( 'f(\theta 1,\theta 2)' );
view(az, el)
