%% Applicazione  di MonteCarlo Approximation e IR ed IS ( importance Sampling Resampling)
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
%
%
% Incrementando il Numero di Particelle e' possibile osservare come
% si riduce l'errore di approssimazione ( teoricamente O(N^-1/2 ))
clear all; close all
mu = 0.8; %
sigma = sqrt(1.5); %
k=1.65; %
c=2; %
num_iterations = 10; %
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
legend('p(x)','q(x)','f(x)*p(x)')


samp_size = 2500;
%SCHEMA IMPORTANCE SAMPLING  eseguito per num_iterations Simulazioni
for t=1:num_iterations
     
    % CAMPIONO da q(x) : X --> Vettore dei samp_size Samples
    Samples = normrnd(mu,sigma,[samp_size 1]);
     
    %Calcolo I Pesi
    j=0;
    for i=1:length(Samples)
        %Considero solo i sample Positivi perche' l'obbiettivo e
        %approssimare l'integrale in [0,5]
        if Samples(i) >=0
            j = j + 1;
            X(j)=Samples(i);
            
        end
    end
    %Calcolo I PESI
    W=p(X )./q(X );
    
    %Calcolo i termini da Sommare nell' approssimazione
    I=W.*f(X);
    %Calcola l'integrale Approssimato: Montecarlo sui campioni Estratti 
    I_hat_is(t) = sum(I)/sum(W);
    
    %     figure(2);
    %     clf;
    %     plot(X,W,'ro');hold on
    
    %% RESAMPLING con Campionamento Inverso Discreto se le particelle 
    %  generate sono poche
    % 
    % controllo se sia necessario ricampionare
    % neff=1/(W*W.');
    % if(neff < samp_size/10)
        
        W=W/sum(W);
        Index=CampionamentoInversoDiscreto(W, samp_size)+1;
        Xres=X(Index);
        %Calcolo I PESI
        Wres=ones(1,samp_size)/samp_size;
        %Calcolo i termini da Sommare nell' approssimazione
        I=Wres.*f(Xres);
        %Calcola l'integrale Approssimato: Montecarlo sui campioni Resampled 
        I_hat_ir(t) = sum(I)/sum(Wres);

        
    %    plot(Xres,Wres,'bo')
    %    legend(['IS-I=' num2str(I_hat_is)],['IR-I=' num2str(I_hat_ir)])
      
    %end
 
    
 
 
end
disp([I_hat_is ; I_hat_ir])
