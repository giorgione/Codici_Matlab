function ParticleEx1
%%
% Particle filter- BootStrap Filter example, adapted from Gordon, Salmond, and Smith paper.
%
% Sistema Dinamico F: R --> R
%
% x = 0.5 * x + 25 * x / (1 + x^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;
%
% Equazione Misura H: R --> R
% y = x^2 / 20 + sqrt(R) * randn;


x = 0.1; % initial state
Q = 1; % process noise covariance
R = 1; % measurement noise covariance
tf = 50; % simulation length

N = 100; % number of particles in the particle filter

xhat = x;
P = 2;
xhatPart = x;

h=figure(1); hold on;
title('Particelle');

%Genero N Particelle al tempo to
% Perturbando casualmente lo stato INIZIALE
for i = 1 : N
    xpart(i) = x + sqrt(P) * randn;
    plot(0,x,'b*');
   
end
plot(zeros(1,N),xpart,'o','Color',[1,0,0]);
  
xArr = [x];
yArr = [x^2 / 20 + sqrt(R) * randn];
xhatArr = [x];
PArr = [P];
xhatPartArr = [xhatPart];

close all;

for k = 1 : tf
    % Evoluzione del Sistema DINAMICO che mi restituisce la predizione
    % dello stato attuale
    x = 0.5 * x + 25 * x / (1 + x^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;
    
    %Misurazioni al tempo attuale: Likelihood
    y = x^2 / 20 + sqrt(R) * randn;
 
    % Extended Kalman filter
    F = 0.5 + 25 * (1 - xhat^2) / (1 + xhat^2)^2;
    P = F * P * F' + Q;
    H = xhat / 10;
    K = P * H' * (H * P * H' + R)^(-1);
    xhat = 0.5 * xhat + 25 * xhat / (1 + xhat^2) + 8 * cos(1.2*(k-1));
    xhat = xhat + K * (y - xhat^2 / 20);
    P = (1 - K * H) * P;
    
    % Particle filter: Genero le Particelle 
    for i = 1 : N
        %Propago le Particelle nel sistema per stimare i PRIOR
        % GENERAZIONE CHE AVVIENE FACENDO EVOLVERE NEL SISTEMA LE
        % PARTICELLE generate al passo precedente
        xpartminus(i) = 0.5 * xpart(i) + 25 * xpart(i) / (1 + xpart(i)^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;
        
        %Misura predetta per la particella secondo il MODELLO DI MISURA:
        %Misura predetta senza la presenza di Rumore
        ypart = xpartminus(i)^2 / 20;
        
        %Errore tra misura Osservata y e misura predetta ypart in funzione della Particella
        vhat = y - ypart;%
        
        %Peso associato alla Particella:  
        %       P( Y(k) ) | P( X(k) ) ) =
        %
        % H(X(k)) + N(0,R) - Yosservato = N(H(X(k))-Yosservato,R)
        %
        %  H(X(k))-Yosservato: ERRORE DI MISURA per la particella(i-sima)
        %  R: matrice di Covarianza dell'errore di misura:    
        %
        % minore è l'errore di misura maggiore è il peso:
        % misura quanto la misura (o osservazione) predetta è VEROSIMILE
        % alla misura precedente
        q(i) = (1 / ( sqrt(R)*sqrt(2*pi)) ) * exp( -vhat^2 / (2*R) );
    end
    
    % Normalize the likelihood of each a priori estimate
    % FACCIO IN MODO CHE i pesi rappresentino una distribuzione (q--> distribuzione)
    qsum = sum(q);    
    q = q / qsum;

    
    colore=q/max(q);
    figure(1);hold on
    for i = 1 : N      
        
        %ridisegno le particelle con i colore 
        
        plot(k,xpart(i),'o','Color',[colore(i) 0 0]);

    end
    % Resample - 
    % RICAMPIONAMENTO delle particelle attraverso un PROCESSO ALEATORIO
    % IDEA:
    % SOPRAVVIVONO le Particelle che hanno PESO MAGGIORE ( errore di misura minimo )
    % VENGONO RIMPIAZZATE le Particelle che hanno PESO MINORE
    
    %Per ogni particella i-esima effettuo il ricampionamento
    for i = 1 : N
        %Genera un numero casuale in [0 1]
        u = rand; % uniform random number between 0 and 1
        qtempsum = 0;
        
        %Somma i pesi fino a quando il totale non è maggiore di u:
        % in questo caso assegno alla  Particella corrente la particella
        % Approssimata
        for j = 1 : N
            
            %Distribuzione cumultativa dei PESI (errore cumulativo sulle misure): cdf(j)
            qtempsum = qtempsum + q(j);
            
            %Nel momento in cui cdf(j>=u) aggiorno la particella xpart(i)
            %con la paricella predetta j-esima
            %
            % L'errore CUMULATIVO ECCEDE LA SOGLIA U
            if qtempsum >= u
                %xpartminus(j): particella predetta che fa eccedere la
                %Soglia di errore u  ---> PARTICELLA CHE HA UN PESO ELEVATO
                xpart(i) = xpartminus(j);
                break;
            end
        end
        
    end
    
    %Disegno le particelle dopo il campionamento
    plot(k*ones(1,N),xpart,'og');
    
    %Calcolo la Media delle particelle dopo il Ricampionamento
    % The particle filter estimate is the mean of the particles. 
    %Lo stato STIMATO è la MEDIA Delle PARTICELLE
    xhatPart = mean(xpart);
    plot(k,xhatPart,'om','MarkerFaceColor','m');
    % Plot the estimated pdf's at a specific time.
    if k == 10
        
        % Particle filter pdf
        pdf = zeros(81,1);
        for m = -40 : 40
            for i = 1 : N
                if  (m <= xpart(i)) && ( xpart(i) < m+1)
                    pdf(m+41) = pdf(m+41) + 1;
                end
            end
        end
        
        figure;
        m = -40 : 40;
        plot(m, pdf / N, 'r');
        hold;
        title('Estimated pdf at k=20');
        disp(['min, max xpart(i) at k = 20: ', num2str(min(xpart)), ', ', num2str(max(xpart))]);
     
        % Kalman filter pdf
        pdf = (1 / sqrt(P) / sqrt(2*pi)) .* exp(-(m - xhat).^2 / 2 / P);
        plot(m, pdf, 'b');
        legend('Particle filter', 'Kalman filter');
        
    end
    
    % Save data in arrays for later plotting
    %Stato
    xArr = [xArr x];
    %Misura
    yArr = [yArr y];
    
    %Stato Stimato da EKF
    xhatArr = [xhatArr xhat];
    
    PArr = [PArr P];    
    %Stato Stimato dal Particle Filter
    xhatPartArr = [xhatPartArr xhatPart];
end

t = 0 : tf;

%figure;
%plot(t, xArr);
%ylabel('true state');

figure;
plot(t, xArr, 'b.', t, xhatArr, 'k-', t, xhatArr-2*sqrt(PArr), 'r:', t, xhatArr+2*sqrt(PArr), 'r:');
axis([0 tf -40 40]);
set(gca,'FontSize',12); set(gcf,'Color','White'); 
xlabel('time step'); ylabel('state');
legend('True state', 'EKF estimate', '95% confidence region'); 

figure;
plot(t, xArr, 'b.', t, xhatPartArr, 'k-');
set(gca,'FontSize',12); set(gcf,'Color','White'); 
xlabel('time step'); ylabel('state');
legend('True state', 'Particle filter estimate'); 

xhatRMS = sqrt((norm(xArr - xhatArr))^2 / tf);
xhatPartRMS = sqrt((norm(xArr - xhatPartArr))^2 / tf);
disp(['Kalman filter RMS error = ', num2str(xhatRMS)]);
disp(['Particle filter RMS error = ', num2str(xhatPartRMS)]);
