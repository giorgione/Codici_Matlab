%Esercizio 1
%
%Tecnica di campionamento Inverso (DISCRETO)

% Simulate the distribution observed in the
% human random digit generation task

% probabilities for each digit
theta = [0.000; ... % digit 0
0.100; ... % digit 1
0.090; ... % digit 2
0.095; ... % digit 3
0.200; ... % digit 4
0.175; ... % digit 5
0.190; ... % digit 6
0.050; ... % digit 7
0.100; ... % digit 8
0.000 ] ... % digit 9
 
% fix the random number generator
seed = 1; rand( 'state' , seed );
 
% let's say we draw K random values
K = 10000;
digitset = 0:9;
Y = randsample(digitset,K,true,theta);
 
 % create a new figure
figure( 1 ); clf;

% Show the histogram of the simulated draws
counts = hist( Y , digitset );
bar( digitset , counts , 'k' );
xlim( [ 0.5 9.5 ] );
xlabel( 'Digit' );
ylabel( 'Frequency' );
title( 'Distribution of simulated draws of human digit generator' );


%Calcolo la Distribuzione cumulativa di theta
thetaC=cumsum(theta);

%genero K campioni Random a distribuzione uniforme in [0 1]
Y1=unifrnd(0,1,1,K);
outcome=zeros(1,K);
%Campionamento Inverso
for i=1:K
    j=1;
    while Y1(i) >= thetaC(j)
        j=j+1;          
    end
    outcome(i)=j-1;
end
% create a new figure
figure(2 ); clf;

% Show the histogram of the simulated draws
counts = hist( outcome , digitset );
bar( digitset , counts , 'k' );
xlim( [ 0.5 9.5 ] );
xlabel( 'Digit' );
ylabel( 'Frequency' );
title( 'Distribution of simulated draws of human digit generator' );

