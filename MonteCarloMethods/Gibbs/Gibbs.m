%Gibbs Sampler da una Gauassiana Bivariata
%% Parameters of the Bivariate normal
u1=0;
u2=0;
mu= [ u1 u2 ];
v=1;
cov=0.7;

sigma = [ v cov; 
          cov v];

%% Initialize the Metropolis sampler
T= 5000; % Set the maximum number of iterations
thetamin = [ -3 -3 ]; % define minimum for theta1 and theta2
thetamax = [ 3 3 ]; % define maximum for theta1 and theta2
seed=1; 
rand( 'state' , seed ); randn('state',seed ); % set the random seed
 
state = zeros( 2 , T ); % Init storage space for the state of the sampler
theta1 = 0; % Start value for theta1
theta2 = 0; % Start value for theta2
t = 1; % initialize iteration at 1
state(1,t) = theta1; % save the current state
state(2,t) = theta2;


%% Start sampling
while t < T % Iterate until we have T samples
    t = t + 1;
    %                    t        t-1   
    %Campiono su P(theta1 | theta2   )
    
    theta1 = normrnd( u1+cov*(theta2-u2),sqrt(1-cov^2) );
    
    
    %                   t        t   
    %Campiono su P(theta2 | theta1   )
    theta2 = normrnd( u2+cov*(theta1-u1),sqrt(1-cov^2) );
    

   

    %% Save state
    state(1,t) = theta1;
    state(2,t) = theta2;
    

end
%% Display histogram of our samples
figure( 1 ); clf;
subplot( 2,2,1 );
nbins = 12;
thetabins1 = linspace( thetamin(1) , thetamax(1) , nbins );
thetabins2 = linspace( thetamin(2) , thetamax(2) , nbins );
hist3( state' , 'Edges' , {thetabins1 thetabins2} );
xlabel( '\theta 1' ); ylabel('\theta 2' ); zlabel( 'counts' );
az = 61; el = 30; view(az, el);

%% Plot the theoretical density
subplot(2,2,2);
nbins = 50;
thetabins1 = linspace( thetamin(1) , thetamax(1) , nbins );
thetabins2 = linspace( thetamin(2) , thetamax(2) , nbins );

[ theta1grid , theta2grid ] = meshgrid( thetabins1 , thetabins2 );
zgrid = mvnpdf( [ theta1grid(:) theta2grid(:)] , mu , sigma );
zgrid = reshape( zgrid , nbins , nbins );
surf( theta1grid , theta2grid , zgrid );
xlabel( '\theta 1' ); ylabel('\theta 2' );
zlabel( 'pdf(\theta 1,\theta 2)' );
view(az, el);
xlim([thetamin(1) thetamax(1)]); ylim([thetamin(2) thetamax(2)]);

subplot(2,2,3);
plot(state(1,1:20),state(2,1:20),'-or')


subplot(2,2,4);
plot(state(1,:),state(2,:),'or')