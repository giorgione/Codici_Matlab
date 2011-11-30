%% Applicazione  di MonteCarlo Approximation e del IS ( importance Sampling )
%  per la stima approssimata del Valore Atteso:
%              I                     1            i  
%  E[ f(x) ] = I   f(x) p(x) dx ~   ---  sum( f( x ) )
%              I                     N
%
%  dove:
%
%  - f(x): funzione arbitraria --> Loss function (in Decision Theoretic)
%
%  - p(x): distribuzione di Probabilita (difficile da campionare)
%
%     i
%  - x : Campioni estratti da p(x)
%
%  Per CAMPIONARE p(x) utilizzo IMPORTANCE SAMPLING :
%
%              ^             1         i      i
%  E[ f(x) ] ~ E [ f(x) ] = ---  sum( w * f( x ) )
%                            N
%  dove:
%
% - g(x): PROPOSAL Distribution  
%         campionata mediante built  --> distribuzione approssimante p(x)
%          in functions MATLAB
%
%    i                                 i      i
% - w : PESI delle PARTICELLE  --> p( x )/g( x )
clear all;
mu = 0.8; %
sigma = sqrt(1.5); %
k=1.65; %
c=2; %
num_iterations = 100; %
% TARGET DISTRIBUTION
p = @(x) x.^(k-1).*exp(-x.^2/2); 

% SAMPLING DISTRIBUTION
q = @(x) c/sqrt(2*pi*sigma^2).*exp(-(mu-x).^2/(2*sigma^2)); 

% FUNZIONE ARBITRARIA sulla quale calcolare il Valore Atteso
f = @(x) 2*sin(pi/1.5*x);

%Considero le distribuzioni nell' intervallo [0 , 5] rispetto al quale
%vogli calcolare l' integrale
fun=@(x)f(x).*p(x);
I_quad = quad(fun,0,5)
 
x=linspace(0,5,100);
figure(1);
plot(x,p(x),'b');hold on
plot(x,q(x),'r');
plot(x,p(x).*f(x),'g'); 
legend('p(x)','q(x)','p(x)*f(x)')

%SCHEMA IMPORTANCE SAMPLING  eseguito per num_iterations Simulazioni
for iter=1:num_iterations
    samp_size = 1000; 
    % CAMPIONO da q(x) : X --> Vettore dei samp_size Samples
    X = normrnd(mu,sigma,[samp_size 1]);
     
    %Calcolo I Pesi
    for i=1:length(X)
        %Considero solo i sample Positivi perche' l'obbiettivo e
        %approssimare l'integrale in [0,5]
        if X(i)>=0
            %Calcolo I PESI
            W(i) = p(X(i))/q(X(i));
            %Calcolo i termini da Sommare nell' approssimazione
            I(i) = W(i)*f(X(i));
        else
            %Assegno peso nullo e termine nell 'integrale nullo in modo da
            %non inficiare l'approssimazione
            W(i)=0;
            I(i)=0;
            %Aggiorna il sample size: numero di samples utilizzati per
            %approssimare l'integrale
            samp_size = samp_size - 1;
        end
    end
    
    %Calcola l'integrale Approssimato (Montecarlo)
    I_hat = sum(I)/sum(W);
    
    [tmp tmp2 non_zero_weights] = find(W);
    %Varianza dei Samples generati
    variance = var(non_zero_weights);
    eff_samp_size = samp_size/(1 + variance);
    %Salvo i seguenti risultati ad ogni Simulazione
    % - Integrale Approssimato
    % - Varianza dei Samples generati
    % - Numero di Campioni Generati
    % - effective sample size
    format short g
    results(iter,:) = [I_hat variance samp_size eff_samp_size];
 
 
end
% Calcolo i Valori MEdi e Varianza dei risultati ottenuti
totals = [mean(results(:,1)) mean(results(:,2)) mean(results(:,3)) mean(results(:,4));...
          var(results(:,1)) var(results(:,2)) var(results(:,3)) var(results(:,4))]