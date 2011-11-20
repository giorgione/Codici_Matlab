function y = cauchy( theta )
%% Returns the unnormalized density of the Cauchy distribution
y = 1 ./ (1 + theta.^2);
