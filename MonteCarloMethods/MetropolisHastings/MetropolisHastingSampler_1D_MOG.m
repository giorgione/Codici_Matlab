%%  Use Metropolis-Hastings procedure to sample from a Mixture of Gaussian
%
%   Utilizzo come PROPOSAL:  
%                  t       t-1              t
%           Q(theta | theta   ) = Gamma(theta| theta*1/tau , tau) 
%
%          
% 
%  Utilizzo starting point differenti per vedere come Campiona la
%  procedura
%
%% Initialize the Metropolis sampler
clc;clear;close all
delete *.eps;
delete *.tex;

T= 500; % Set the maximum number of iterations
tau = 0.25:0.25:9; % Set standard deviation of normal proposal density
tau=tau.';

thetamin = 1; 
thetamax = 100; % define a range for starting values

Mixture=@(X,m1,m2,v1,v2,p)...
    p(1)*exp(-(X-m1).^2./(2*v1.^2))./(sqrt(2*pi)*v1)+...
    p(2)*exp(-(X-m2).^2./(2*v2.^2))./(sqrt(2*pi)*v2);

P=@(X)Mixture(X,30,50,4,3,[0.3,0.7]);

%Genero K processi di Campionamento
K=length(tau);
theta = zeros( K , T ); % Init storage space for our samples
Accept = zeros( K , T ); % Statistics for Acceptance

seed=1; rand( 'state' , seed ); randn('state',seed ); % set the random seed
%Genero il punto Iniziale
theta(:,1) = unifrnd( thetamin , thetamax,K,1); % Generate start value


%% Start sampling
t = 1;
while t < T % Iterate until we have T samples
    t = t + 1;
    
    % Propose a new value for theta using a GAMMA PROPOSAL DENSITY
    thetastar = gamrnd( theta(:,t-1).*tau ,ones(K,1)./tau);
    
    %Calcolo i fattori che caratterizzano la ACCCEPTANCE PROBABILITY
    alpha1=gampdf(theta(:,t-1),thetastar.*tau ,ones(K,1)./tau)./gampdf(thetastar,theta(:,t-1).*tau ,ones(K,1)./tau);
    
    
    
    % Calcolo del ACCEPTANCE RATIO:
    %                  *          t-1     *
    %          P( theta  ) Q(theta | theta)
    % min([1,---------------------------])
    %                t-1          *       t-1
    %         P( theta  )  Q(theta | theta   )
    
    alpha0=P(thetastar)./ P(theta(:,t-1));

    % Calculate the ACCEPTANCE RATIO
    alpha = min( [ ones(K,1)  alpha0.*alpha1].' );
    % Draw a UNIFORM DEVIATE from [ 0 1 ]
    u = rand(1,K);
    
    for i=1:K
        % Do we accept this proposal?
        if u(i) < alpha(i)
            theta(i,t) = thetastar(i); % If so, proposal becomes new state
            Accept(i,t)=1;
        else
            theta(i,t) = theta(i,t-1); % If not, copy old state
            Accept(i,t)=0;
        end
    end
end

nbins = 200;
thetabins = linspace( thetamin , thetamax , nbins );
%% Display histogram of our samples
for i=1:K
    
    figure( i ); clf;
   
    subplot( 3,1,1 );


    counts = hist( theta(i,:) , thetabins );
    bar( thetabins , counts/sum(counts) , 'g' );
    xlim( [ thetamin thetamax ] );
    xlabel( '$\theta$','Interpreter', 'latex'); 
    ylabel( '$p(\theta)$','Interpreter', 'latex');
    
    %% Overlay the theoretical density
    y = P( thetabins);
    hold on;
    plot( thetabins , y/sum(y) , 'r-' , 'LineWidth' , 3 );
    set( gca , 'YTick' , [] );
    %legend('Distribuzione Campioni Estratti','Distr. Campionata')

    %% Display history of our samples
    subplot( 3,1,2:3 );
    stairs(1:T, theta(i,:) ,  'g-' );
    xlabel( 't'); ylabel( '$\theta^t$','Interpreter', 'latex');
    %legend('Processo Sampling nel Tempo')
    ylim( [ thetamin thetamax ] );
    %set( gca , 'YDir' , 'reverse' );
    
    fileimg=['MetropolisHastingExample' num2str(i) '.eps'];
    
    saveas(gcf,fileimg, 'psc2')
    close gcf
end

%Sovrappongo tutte i processi di Campionamento
Tvec=1:T;
Tvec=repmat(Tvec',1,K);
figure(K+1);clf
plot( Tvec,theta.'  , '-');
xlabel( 't' ); ylabel( '\theta' );
set( gca , 'YDir' , 'reverse' );
%legend(['Sampling t_0=' num2str(theta(1,1))],...
%       ['Sampling t_0=' num2str(theta(2,1))],...
%       ['Sampling t_0=' num2str(theta(3,1))],...
%       ['Sampling t_0=' num2str(theta(4,1))]);
%set( gca , 'YDir' , 'reverse' );

A=sum(Accept,2);

FID = fopen('MetropolisHastingSimulation.tex', 'w');

fprintf(FID, '\\begin{center}\n');
fprintf(FID, '\\begin{table}\\label{tb: Simulation Metropolis}\n');
fprintf(FID, '\\begin{tabular}{|ccccc|}\\hline \n');
fprintf(FID, 'Simulation & $\\theta_0$ & $\\sigma$ & Acceptance ratio & Reject ratio)\\\\ \\hline \n');

%Generate some statistic
for i=1:K
    display('---------------------------------------');
    display(['PROCESS ' num2str(i)]);
    display(['To: ' num2str(theta(i,1))]);
    display(['sigma: ' num2str(tau(i))]);
    display(['Total number of iteration: ' num2str(T)]);
    display(['Accepted Samples: ' num2str(A(i)/T)]);
    display(['Rejected Samples: ' num2str(1-A(i)/T)]);
    display('---------------------------------------');
    legenda{i}=['\theta_0=' num2str(theta(1,1)) 'sigma=' num2str(tau(i))];
  
    fprintf(FID, '%d & %8.2f & %8.2f & %8.2f & %8.2f \\\\  \hline',i, theta(i,1), tau(i),A(i)/T , 1-A(i)/T);
    if i==length(K)
        fprintf(FID, '\\hline ');
    end
    fprintf(FID, '\n');
    
end
fprintf(FID, '\\end{tabular}\n\\caption{ Metropolis simulation results varying $\\theta_0$ and $\\sigma$ }\n');
fprintf(FID, '\\end{table}\n');
fprintf(FID, '\\end{center}\n');

fclose(FID);
legend(legenda)


    