%% Chapter 2. Metropolis sampler for Mallows model
% samples orderings from a distribution over orderings

%% Initialize model parameters
lambda = 0.1; % scaling parameter
labels = { 'Washington' , 'Adams' , 'Jefferson' , 'Madison' , 'Monroe' };
omega  = [ 1 2 3 4 5 ]; % correct ordering
L      = length( omega ); % number of items in ordering

%% Initialize the Metropolis sampler
T     = 500; % Set the maximum number of iterations
seed=1; rand( 'state' , seed ); randn('state',seed ); % set the random seed
theta = zeros( L , T ); % Init storage space for our samples
theta(:,1) = randperm( L ); % Random ordering to start with

%% Start sampling
t = 1; 
while t < T % Iterate until we have T samples
    t = t + 1;       
    
    % Our proposal is the last ordering but with two items switched
    lasttheta = theta(:,t-1); % Get the last theta  
    % Propose two items to switch
    whswap = randperm( L ); whswap = whswap(1:2);      
    theta_star = lasttheta;
    theta_star( whswap(1)) = lasttheta( whswap(2));
    theta_star( whswap(2)) = lasttheta( whswap(1));
     
    % calculate Kendall tau distances
    dist1 = kendalltau( theta_star , omega );
    dist2 = kendalltau( lasttheta  , omega );
    
    % Calculate the acceptance ratio 
    pratio = exp( -dist1*lambda ) / exp( -dist2*lambda );
    alpha = min( [ 1 pratio ] );        
    u = rand; % Draw a uniform deviate from [ 0 1 ]    
    if u < alpha  % Do we accept this proposal?
        theta(:,t) = theta_star; % proposal becomes new theta       
    else
        theta(:,t) = lasttheta; % copy old theta
    end
    % Occasionally print the current state
    if mod( t,10 ) == 0
        fprintf( 't=%3d\t' , t );
        for j=1:L
            fprintf( '%15s' , labels{ theta(j,t)} );
        end
        fprintf( '\n' );
    end
end

