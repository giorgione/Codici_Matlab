%% Chapter 2. Exercise 5
%  Use Metropolis-Hastings procedure to sample from WEIBULL DENSITY
%  Utilizzo starting point differenti per vedere come Campiona la
%  procedura
%
%% Initialize the Metropolis sampler
T= 500; % Set the maximum number of iterations
sigma = 3; % Set standard deviation of normal proposal density
 
thetamin = -30; 
thetamax = 30; % define a range for starting values

%Genero K processi di Campionamento
K=4;
theta = zeros( K , T); % Init storage space for our samples
seed=1; rand( 'state' , seed ); randn('state',seed ); % set the random seed
%Genero il punto Iniziale
theta(:,1) = round(unifrnd( thetamin , thetamax,K,1)); % Generate start value


%% Start sampling
%parametro della gamma
tau=1;
%Parametri della weibull distribution
a=2;
b=1.9;
t = 1;
Acceptance=zeros(1,K);
while t < T % Iterate until we have T samples
    t = t + 1;
    
    % Propose a new value for theta using a GAMMA PROPOSAL DENSITY
    thetastar = gamrnd( theta(:,t-1)*tau ,ones(K,1)./tau);
    
    %Calcolo i fattori ce caratterizzano la ACCCEPTANCE PROBABILITY
    alpha1=gampdf(theta(:,t-1),thetastar*tau ,1/tau)./gampdf(thetastar,theta(:,t-1)*tau ,1/tau);
    alpha0=wblpdf(thetastar,a,b)./wblpdf(theta(:,t-1),a,b);
    
    % Calculate the acceptance ratio
    alpha = min( [ ones(K,1)  alpha0.*alpha1].' );
                                
                            
    % Draw a uniform deviate from [ 0 1 ]
    u = rand(1,K);
    
    for i=1:K
        % Do we accept this proposal?
        if u(i) < alpha(i)
            theta(i,t) = thetastar(i); % If so, proposal becomes new state
            Acceptance(i)=Acceptance(i)+1;
        else
            theta(i,t) = theta(i,t-1); % If not, copy old state
        end
    end
end
Acceptance=Acceptance/T;

nbins = 200;
thetabins = linspace( thetamin , thetamax , nbins );
%% Display histogram of our samples
for i=1:K
    
    figure( i ); clf;
    subplot( 3,1,1 );


    counts = hist( theta(i,:) , thetabins );
    bar( thetabins , counts/sum(counts) , 'k' );
    xlim( [ thetamin thetamax ] );
    xlabel( '\theta' ); ylabel( 'p(\theta)' );

    %% Overlay the theoretical density
    y = wblpdf( thetabins,a,b );
    hold on;
    plot( thetabins , y/sum(y) , 'r-' , 'LineWidth' , 3 );
    set( gca , 'YTick' , [] );

    %% Display history of our samples
    subplot( 3,1,2:3 );
    stairs( theta(i,:) , 1:T , 'k-' );
    ylabel( 't' ); xlabel( '\theta' );
    xlim( [ thetamin thetamax ] );
    set( gca , 'YDir' , 'reverse' );
end

%Sovrappongo tutte i processi di Campionamento
T=1:T;
T=repmat(T',1,K);
figure(K+1);clf
plot( theta.' , T , '-');
ylabel( 't' ); xlabel( '\theta' );
set( gca , 'YDir' , 'reverse' );