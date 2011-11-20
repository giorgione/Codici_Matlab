%% Chapter 2. Metropolis procedure to sample from Bivariate Exponential
% BLOCKWISE UPDATING. Use a uniform proposal distribution
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

%% Start sampling
t = 1;
while t < T % Iterate until we have T samples
    t = t + 1;
    % Propose a new value for theta
    thetastar = unifrnd( thetamin , thetamax );
    pratio = bivexp( thetastar(1) , thetastar(2) ) / ...
    bivexp( theta(1,t-1), theta(2,t-1) );
    alpha = min( [ 1 pratio ] ); % Calculate the acceptance ratio
    u = rand; % Draw a uniform deviate from [ 0 1 ]
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
