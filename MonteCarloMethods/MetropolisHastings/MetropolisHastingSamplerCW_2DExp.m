%% Chapter 2.
%Metropolis procedure to sample from Bivariate Exponential
% Component−wise updating. Use a normal proposal distribution
clc;clear;close all
%% Initialize the Metropolis sampler
T= 5000; % Set the maximum number of iterations
thetamin = [ 0 0 ]; % define minimum for theta1 and theta2
thetamax = [ 8 8 ]; % define maximum for theta1 and theta2
seed=1; 
rand( 'state' , seed ); randn('state',seed ); % set the random seed
 
state = zeros( 2 , T ); % Init storage space for the state of the sampler
theta1 = unifrnd( thetamin(1) , thetamax(1), 1,1 ); % Start value for theta1
theta2 = unifrnd( thetamin(2) , thetamax(2),1,1 ); % Start value for theta2

t = 1; % initialize iteration at 1
state(1,t) = theta1; % save the current state
state(2,t) = theta2;
Acceptance=zeros(1,2);
%% Start sampling
while t < T % Iterate until we have T samples
    t = t + 1;

    %% Propose a new value for theta1
    newtheta1 = unifrnd( thetamin(1) , thetamax(1),1,1 );
    pratio = bivexp( newtheta1,theta2)/bivexp( theta1,theta2);
    
    alpha = min( [ 1 pratio ] ); % Calculate the acceptance ratio
    u = rand; % Draw a uniform deviate from [ 0 1 ]
    if u < alpha % Do we accept this proposal?
        theta1 = newtheta1; % proposal becomes new value for theta1
        Acceptance(1)=Acceptance(1)+1;
    end

    %% Propose a new value for theta2
    newtheta2 = unifrnd( thetamin(2) , thetamax(2),1,1 );
    pratio = bivexp( newtheta2,theta1)/bivexp( theta2,theta1);
    alpha = min( [ 1 pratio ] ); % Calculate the acceptance ratio
    u = rand; % Draw a uniform deviate from [ 0 1 ]
    if u < alpha % Do we accept this proposal?
        theta2 = newtheta2; % proposal becomes new value for theta2
        Acceptance(2)=Acceptance(2)+1;
    end

    %% Save state
    state(1,t) = theta1;
    state(2,t) = theta2;

end

%% Display histogram of our samples
figure( 1 ); clf;
subplot( 1,2,1 );
nbins = 12;
thetabins1 = linspace( thetamin(1) , thetamax(1) , nbins );
thetabins2 = linspace( thetamin(2) , thetamax(2) , nbins );
hist3( state' , 'Edges' , {thetabins1 thetabins2} );
xlabel( '\theta 1' ); ylabel('\theta 2' ); zlabel( 'counts' );
az = 61; el = 30; view(az, el);

%% Plot the theoretical density
subplot(1,2,2);
nbins = 50;
thetabins1 = linspace( thetamin(1) , thetamax(1) , nbins );
thetabins2 = linspace( thetamin(2) , thetamax(2) , nbins );
[ theta1grid , theta2grid ] = meshgrid( thetabins1 , thetabins2 );
zgrid = bivexp(  theta1grid(:) ,theta2grid(:));
zgrid = reshape( zgrid , nbins , nbins );
surf( theta1grid , theta2grid , zgrid );
xlabel( '\theta 1' ); ylabel('\theta 2' );
zlabel( 'pdf(\theta 1,\theta 2)' );
view(az, el);
xlim([thetamin(1) thetamax(1)]); ylim([thetamin(2) thetamax(2)]);
