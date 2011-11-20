% Chapter 2. Use Metropolis procedure to sample from Cauchy density

%% Initialize the Metropolis sampler
T= 500; % Set the maximum number of iterations
sigma = 1; % Set standard deviation of normal proposal density
thetamin = -30; 
thetamax = 30; % define a range for starting values
theta = zeros( 1 , T ); % Init storage space for our samples
seed=1; rand( 'state' , seed ); randn('state',seed ); % set the random seed
theta(1) = unifrnd( thetamin , thetamax ); % Generate start value


%% Start sampling
t = 1;
while t < T % Iterate until we have T samples
    t = t + 1;
    % Propose a new value for theta using a NORMAL PROPOSAL DENSITY
    thetastar = normrnd( theta(t-1) , sigma );
    
    % Calculate the acceptance ratio
    alpha = min( [ 1 cauchy( thetastar ) / cauchy( theta(t-1) ) ] );
    % Draw a uniform deviate from [ 0 1 ]
    u = rand;
    
    % Do we accept this proposal?
    if u < alpha
        theta(t) = thetastar; % If so, proposal becomes new state
    else
        theta(t) = theta(t-1); % If not, copy old state
    end
end


%% Display histogram of our samples
figure( 1 ); clf;
subplot( 3,1,1 );
nbins = 200;
thetabins = linspace( thetamin , thetamax , nbins );
counts = hist( theta , thetabins );
bar( thetabins , counts/sum(counts) , 'k' );
xlim( [ thetamin thetamax ] );
xlabel( '\theta' ); ylabel( 'p(\theta)' );

%% Overlay the theoretical density
y = cauchy( thetabins );
hold on;
plot( thetabins , y/sum(y) , 'r-' , 'LineWidth' , 3 );
set( gca , 'YTick' , [] );

%% Display history of our samples
subplot( 3,1,2:3 );
stairs( theta , 1:T , 'k-' );
ylabel( 't' ); xlabel( '\theta' );
set( gca , 'YDir' , 'reverse' );
