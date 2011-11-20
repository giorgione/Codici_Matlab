%% Explore the Normal distribution N( mu , sigma )
mu    = 100; % the mean
sigma = 15;  % the standard deviation
xmin  = 70;  % minimum x value for pdf and cdf plot
xmax  = 130; % maximum x value for pdf and cdf plot
n     = 100; % number of points on pdf and cdf plot
k     = 10000; % number of random draws for histogram

% create a set of values ranging from xmin to xmax
x = linspace( xmin , xmax , n ); 
p = normpdf( x , mu , sigma ); % calculate the pdf 
c = normcdf( x , mu , sigma ); % calculate the cdf 

figure( 1 ); clf; % create a new figure and clear the contents

subplot( 1,3,1 );
plot( x , p , 'k-' );
xlabel( 'x' ); ylabel( 'pdf' );
title( 'Probability Density Function' );

subplot( 1,3,2 );
plot( x , c , 'k-' );
xlabel( 'x' ); ylabel( 'cdf' );
title( 'Cumulative Density Function' );

% draw k random numbers from a N( mu , sigma ) distribution
y = normrnd( mu , sigma , k , 1 ); 

subplot( 1,3,3 );
hist( y , 20 );
xlabel( 'x' ); ylabel( 'frequency' );
title( 'Histogram of random values' );