function y = bivexp( theta1 , theta2 )
%% Returns the density of a bivariate exponential function
lambda1 = 0.5; % Set up some constants
lambda2 = 0.1;
lambda = 0.01;
maxval = 8;
y = exp( -(lambda1+lambda)*theta1-(lambda2+lambda)*theta2-lambda*maxval );
