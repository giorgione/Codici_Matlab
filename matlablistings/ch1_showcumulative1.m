%% Create illustration for the cumulative distribution

% probabilities for each digit
theta = [0.000;  ... % digit 0
         0.100;  ... % digit 1
         0.090;  ... % digit 2
         0.095;  ... % digit 3
         0.200;  ... % digit 4
         0.175;  ... % digit 5
         0.190;  ... % digit 6
         0.050;  ... % digit 7
         0.100;  ... % digit 8
         0.000 ] ...  % digit 9 

c = cumsum( theta );
figure( 1 ); clf;

digitset = 0:9;
stairs( digitset , c , 'k' );
xlim( [ 0 9 ] );
ylim( [ 0 1 ] );
ylabel( 'F(X)' );
xlabel( 'X' );

hold on;
U = 0.8;
whdigit = digitset( find( c > U , 1 , 'first' ));
%plot( [ 0 whdigit ] , [ U U ] , 'k--' );
%plot( [ whdigit whdigit ] , [ U 0 ] , 'k--' );



