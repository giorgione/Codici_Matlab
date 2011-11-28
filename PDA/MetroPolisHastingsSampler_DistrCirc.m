%function [Samples]=MetroPolisHastingsSampler(P,Theta,Q)
%
% Campionamento MetropolisHastings per una Distribuzione Multivariata.
% Parametri:
%   
%   P: Distribuzione che desidero campionare come funzione Vettoriale dei
%      parametri
%
%   Theta: Vettore delle varibili su cui condizionare P
%                                                                       t-1
%   Q: Proposal Distribution come funzione del campione precedente Theta
%              t       t-1
%       Q(theta | Theta   )
%
%   N: Numero Campioni
function [Theta]=MetroPolisHastingsSampler_DistrCirc(P,Theta,NTheta,T,Sigma)

 
%Fisso la proposal come una Bivariate Normal Distribution di parametri a e b
Q=@(x,a,b)mvnpdf(x,a,b);
Qsample=@(a,b)mvnrnd(a,b);
 
 
ThetaSet=1:NTheta;

for t=2:T
   
    %Genera i Campioni per ogni singola Variabile presente nella
    %Distribuzione Congiunta
    for i=1:NTheta
        ThetaNot=setdiff(ThetaSet,i);
       
        ThetaOld=Theta(i,t-1);
        %                                     t
        %Genero il Nuovo Campione per Theta(i)
        ThetaStar=Qsample(ThetaOld,Sigma);
        
        %Calcolo il rapporto 
        %     t-1    *          *      t
        % q( a    | a  ) /  q( a    | a  ) 
        Proposal_ratio = Q( ThetaStar  , ThetaOld  , Sigma)/Q( Theta(i,t-1)  , ThetaStar  ,Sigma );
      
        %Calcolo il rapporto 
        %     *      t-1          t-1    t-1
        % p( a    , b    ) /  p( a    | b    ) 
        %
        % rapporto delle Likelihood * rapporto dei prior
        ThetaParams=zeros(NTheta,1);
        ThetaParams(i)=ThetaStar;
        ThetaParams(ThetaNot)=Theta(ThetaNot,t-1);
        P_ratio_num = P( ThetaParams);        
        ThetaParams(i)=ThetaOld;
        P_ratio_denum =P( ThetaParams);       
        
        P_ratio=P_ratio_num/P_ratio_denum;
        
        % Accept or reject 
        alpha = min( [ 1 P_ratio*Proposal_ratio ] );   
        u = rand; % Draw a uniform deviate from [ 0 1 ]  
        if u < alpha % Do we accept this proposal?
            Theta(i,t)=ThetaStar; % If so, proposal becomes new state for param A
            accept( 1,t ) = 1;
        else
            Theta(i,t)=Theta(i,t-1);

        end
       
    end
    
end
Samples=Theta;
end
