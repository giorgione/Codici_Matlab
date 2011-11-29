%% Chapter 2. 
% Utilizzo Metropolis Sampler per campionare la Cauchy density
% caratterizzata da:
% 
%                1
% P(theta)= --------------
%           pi*(1+theta^2)      
%
%                                                
% Utilizzo come PROPOSAL:  
%                  t       t-1                  t
%           Q(theta | theta   ) = N(theta| theta , var) 
%
%           Gaussiana Centrata nella Stima attuale di theta
%
%  Utilizzo starting point differenti per vedere com lavora l'algoritmo in
%  funzione dello STARTING POINT INIZIALE
%  procedura
%
%% Initialize the Metropolis sampler
T= 500; % Set the maximum number of iterations
sigma = 1; % Set standard deviation of normal proposal density
sigmao= 10;
thetamin = -30; 
thetamax = 30; % define a range for starting values

%Genero K processi di Campionamento
K=4;
theta = zeros( K , T ); % Init storage space for our samples
seed=1; rand( 'state' , seed ); randn('state',seed ); % set the random seed
%Genero il punto Iniziale
theta(:,1) = unifrnd( thetamin , thetamax,K,1); % Generate start value


%% Start sampling
t = 1;
while t < T % Iterate until we have T samples
    t = t + 1;
    
    %% Propose a new value for theta using a NORMAL PROPOSAL DENSITY Q
    thetastar = normrnd( theta(t-1) , sigma,K,1);
    
    %% Calcolo del ACCEPTANCE RATIO:
    %                  t
    %          P( theta  )
    % min([1,--------------])
    %                t-1
    %         P( theta  )
    alpha = min( [ ones(K,1) cauchy( thetastar )./ cauchy( theta(:,t-1) ) ].' );
    
    %% Draw a UNIFORM DEVIATE from [ 0 1 ]
    u = rand(1,K);
    
    for i=1:K
        % Do we accept this proposal?
        if u(i) < alpha(i)
            theta(i,t) = thetastar(i); % If so, proposal becomes new state
        else
            theta(i,t) = theta(i,t-1); % If not, copy old state
        end
    end
end

tnbins = 200;
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
    y = cauchy( thetabins );
    hold on;
    plot( thetabins , y/sum(y) , 'r-' , 'LineWidth' , 3 );
    set( gca , 'YTick' , [] );

    %% Display history of our samples
    subplot( 3,1,2:3 );
    stairs( theta(i,:) , 1:T , 'k-' );
    ylabel( 't' ); xlabel( '\theta' );
    set( gca , 'YDir' , 'reverse' );
end

%Sovrappongo tutte i processi di Campionamento
T=1:T;
T=repmat(T',1,K);
figure(K+1);clf
plot( theta.' , T , '-');
ylabel( 't' ); xlabel( '\theta' );
set( gca , 'YDir' , 'reverse' );
