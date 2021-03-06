%function [Z, corres, corres_car] = MCMCSamplesJointStatesParametrerizationWCar...
%    (prevZ, X, Xc, Y, KLT, corres, corres_car, params, fn)
%
% MCMC  Sampling Procedure
%
% INPUT
% - prevZ: Estimated Status of Variables in the Model from the previous
%          iteration
% - X: PERSON observation
% - Xc: CAR observation
% - Y:
% - KLT: GROUND FEATURES
% - corres:  
% - corres_car: Corrispondences for Car
% - params:  Model Parameters
% - fn :   Current Frame Index
%
% OUTPUT: 
%   - Z
%   - corres:
%   - corres_car:
%
% 
function [Z, corres, corres_car] = MCMCSamplesJointStatesParametrerizationWCar...
    (prevZ, X, Xc, Y, KLT, corres, corres_car, params, fn)

%% Initialization of VARIABLES PARAMTERS: COVARIANCE MATRIX-PRECISION
[nsamples, N, Y, Z, prevZ, newdet, newcdet, newgfeats, corres, corres_car] =...
    initVariables(prevZ, X, Y, Xc, KLT, corres, corres_car, params);

%Contiene Statistiche sui sample generati per ogni tid variable
samplesrecord = zeros(1, prevZ.nCarTargets + length(newcdet) + prevZ.nTargets + 1 + length(newdet) + prevZ.nFeats);
%Contiene Statistiche sui samples accettati per ogni tid variable
acceptrecord = zeros(1, prevZ.nCarTargets + length(newcdet) + prevZ.nTargets + 1 + length(newdet) + prevZ.nFeats);

%Memorizzo i Valori delle Probabilita'
probrecord1 = [1]; 
probrecord2 = [];

nantrial = [];
% initialize the sample
% sample from prevZ
% maybe we can go for multiple times to guarantee multi-modality
TrialFig=1:params.nretry;
TrialFig=20+TrialFig;
for trial = 1:params.nretry
    tempcnt = 0;     
    figure(TrialFig(trial)); 
    title(['Trial ' num2str(trial)]);grid on; hold on;
    Colore=rand(1,3);
    
    %Genero una Variable a caso nel Modello
    initidx = ceil(rand * prevZ.nSamples);
    
    %Get the FIRST SAMPLE: STARTING POINT
    [sample] = getInitialSample(prevZ, Z, X, Xc, KLT, newdet, newcdet, newgfeats, params, initidx);    
    %DrawSample(sample,TrialFig(trial),ColorTrial) 
    % Disegno la posizione della Camera
    plot3(sample.cam(7,:),sample.cam(8,:),sample.cam(2,:),'o','MarkerFaceColor',Colore,'MarkerEdgeColor',Colore,'MarkerSize',5);    %Get Priors 
    % Disegno la posizione dei TARGET
    plot3(sample.car(1,:),sample.car(2,:),sample.car(3,:),'d','MarkerFaceColor',Colore,'MarkerEdgeColor',Colore,'MarkerSize',5);    %Get Priors 
     
    plot3(sample.per(1,:),sample.per(2,:),sample.per(3,:),'s','MarkerFaceColor',Colore,'MarkerEdgeColor',Colore,'MarkerSize',5);    %Get Priors 
    
    [prior, prior_car] = initPriors(sample, prevZ, newdet, newcdet, params);
    if (sum(sum(isnan(prior))) > 0) || (sum(sum(prior == inf)) > 0) || (sum(sum(prior <= 0)) > 0)
        keyboard;
    end
    
    %%%% UPDATE COVARIANCE MATRIX for Each VARIABLE in the Model at current time:
    %
    % - params.Pert{1} = CAMERA COV. MATRIX
    %
    % - params.Pert{[2:prevZ.nTargets][prevZ.nTargets:prevZ.nTargets+length(newdet)]}
    %                  old detections             new detections 
    %
    % - params.cPert{[1:prevZ.nCarTargets][prevZ.nCarTargets:length(newcdet)]}
    %                  old detections         new detections
    %
    [params] = getProposalCovariance(prevZ, initidx, newdet, newcdet, params);
    
    %% Variables to be Sampled
    tempcam = [];    
    tempper = [];    
    tempcar = [];    
    tempbeta = [];   
    tempgfeat = [];
    
    
    if isfield(params, 'camscnt')
        % For each Target how many sample do we need
        %                         
        perTargSample = [params.camscnt * params.thinning, ... %camera
                 params.perscnt * params.thinning * ones(1, prevZ.nTargets + length(newdet)),... %PERSON
                 params.featscnt * params.thinning * ones(1, prevZ.nFeats),.... %GROUND FEAT
                 params.carscnt * params.thinning * ones(1, prevZ.nCarTargets + length(newcdet))]; %CAR
        
             %Cumulative Distribution
        cumSamples = cumsum(perTargSample) / sum(perTargSample);
        N = params.burnin + sum(perTargSample);
    end
    
    %Calcolo la LOG POSTERIOR per lo STARTING POINT
    probrecord1 = [];
    probrecord1 = computeLogPosterior(sample, prevZ, X, Y, Xc, KLT, corres, corres_car, params);
    
    max_prob = -Inf;
    maxsample = sample;    
    temptempsamples = [];
    
    %Total Number of ITERATION of SAMPLING
    for i = 1:N
        dbgcnt = 0;
        tempcnt = tempcnt + 1;
      
        % SELECT  which Sampling VARIABLE (CAMERA-PERSON-CAR-GROUND):
        %
        % tid: index of target variable to be sampled
        % Devo generare campioni per ogni id nell' insieme 
        %
        %          [1,2,...,NVariables] 
        %dove:
        % { 1 ,    --> Variabile Camera 
        %   2
        %   ..     --> Targets Persone gia' presenti nel modello
        %   ..          
        %   2+prevZ.nTargets,                
        %   2+prevZ.nTargets+1
        %   ..                  --> Nuove Persone rilevate dal detector
        %   ..
        %   2+prevZ.nTargets+1+length(newdet), 
        %   2+prevZ.nTargets+length(newdet)+1,
        %   ..                   --->  Auto presenti nel Modello
        %   ..
        %   2+prevZ.nTargets+length(newdet)+1+prevZ.nCarTargets 
        %   2+prevZ.nTargets+length(newdet)+1+prevZ.nCarTargets+1 
        %   ..
        %   ..               --->  Nuove Auto detettate
        %   2+prevZ.nTargets+length(newdet)+1+prevZ.nCarTargets+1 
        % 
  
        if isfield(params, 'camscnt')
            %Number of Total Variables to be Sampled
            NVariables=prevZ.nFeats + prevZ.nTargets + length(newdet) + prevZ.nCarTargets + length(newcdet) + 1;
          
            tid = ceil(rand * NVariables);
            tid = sum(rand > cumSamples) + 1;
        else
            %escludi auto 
            tid = ceil(rand * (prevZ.nFeats + prevZ.nTargets + length(newdet) + 1));
        end
        
        oid = 0; 
        sid = 0;
        %tsample e' settato con il campione precedente
        tsample = sample;
        
        
        if tid == 1
            %CAMERA SAMPLING
            % Object ID = 0
            oid = 0; 
            % Variable ID
            sid = tid;
            tsample.cam = tsample.cam + mvnrnd(zeros(1, params.ncamv), params.Pert{1})';
            if (params.useCamAnnealing == 1)
                params.Pert{1} = params.Pert{1} .* params.AnnealingConst;
            end
            
        elseif tid <= 1+size(sample.per, 2)  %PERSON SAMPLING
            % Object ID = 1
            oid = 1;
            %Id of Variable to be Sampled
            sid = tid - 1; 
            
            %Current value for the sample
            temp = tsample.per(end, tid-1);
            if rand < params.flipObj
                temp = ~temp;
            end
            %New sample for the Variable tid-1
            tsample.per(:, tid-1) =...
            [tsample.per(1:params.nperv, tid-1) + mvnrnd(zeros(1, params.nperv), params.Pert{tid})';
                temp];
            
            % BETA SAMLING: intercation between PERSONs already present in
            % the Model and Variable to be Sampled (tid-1)
            for j = 1:(tid-2)
                temp = rand;
                if (j == tid - 1) || temp < 0.8, continue; end;
                if temp > 0.9
                    tsample.beta(j, tid - 1, 1) = 1;
                    tsample.beta(j, tid - 1, 2) = 0;
                else
                    tsample.beta(j, tid - 1, 1) = 0;
                    tsample.beta(j, tid - 1, 2) = 1;
                end
            end
            % intercation between Variable to be Sampled (tid-1) 
            % and New detection - PERSONs 
            for j = tid:(prevZ.nTargets + length(newdet))
                temp = rand;
                if (j == tid - 1) || temp < 0.8, continue; end;
               
                if temp > 0.9
                    tsample.beta(tid - 1, j, 1) = 1;
                    tsample.beta(tid - 1, j, 2) = 0;
                else
                    tsample.beta(tid - 1, j, 1) = 0;
                    tsample.beta(tid - 1, j, 2) = 1;
                end
            end
            
        elseif tid <= 1 + size(sample.per, 2) + prevZ.nFeats  %GROUND FEATURE SAMPLING
            % Object ID = 2
            oid = 2; 
            sid = tid - 1 - size(sample.per, 2);
            
            temp = tsample.gfeat(end, sid);
            if rand < params.flipFeat
                temp = ~temp;
            end
            tsample.gfeat(:, sid) = [tsample.gfeat(1:(params.ngfeat-1), sid) + mvnrnd(zeros(1, params.ngfeat-1), params.PertGF{sid})'; temp];
        else %CAR
             % Object ID = 3
            oid = 3; 
            sid = tid - 1 - size(sample.per, 2) - prevZ.nFeats;
            
            temp = tsample.car(end, sid);
            if rand < params.flipObj
                temp = ~temp;
            end
            tsample.car(:, sid) = [tsample.car(1:params.ncarv, sid) + mvnrnd(zeros(1, params.ncarv), params.cPert{sid})'; temp];
        end
        
        % evaluate acceptance ratio
        [ar, tprior, tprior_car] = getAcceptanceRatio(prevZ, X, Y, Xc, corres, corres_car, KLT, sample, tsample, oid, sid, params, prior, prior_car); % compute acceptance ratio
        if isnan(ar)
%             keyboard
        end

        % accepted sample tsample: Contiene il campione Congiunto delle
        % VARIABILI 
        if ar > rand
            %Store the current temporary sample
            sample = tsample;  
         
            
            prior = tprior;
            prior_car = tprior_car;
            
            %Increment the Count for the Accpeted Samples relative to the tid Variable
            acceptrecord(tid) = acceptrecord(tid) + 1;
           
            %Store the Cumulative Probability
            if length(probrecord1) == 0
                probrecord1 = log(ar);
            else
                probrecord1 = [probrecord1, probrecord1(end) + log(ar)];
            end
            
            if probrecord1(end) > max_prob
                maxsample = sample;
                max_prob = probrecord1(end);
            end
        end
        
        %Increment the Count for the Generated Samples relative to the tid Variable
        samplesrecord(tid) = samplesrecord(tid) + 1;
        if i > params.burnin
            if mod(i - params.burnin, params.thinning) ~= 0, continue; end
            % save them.
            tempcam(:, end + 1) = getNextCamera(sample.cam, params);
          
            tempper(:, end + 1) = reshape([params.Aper * sample.per(1:params.nperv, :); sample.per(params.nperv + 1, :)], prod(size(sample.per)), 1);
            tempcar(:, end + 1) = reshape([params.Acar * sample.car(1:params.ncarv, :); sample.car(params.ncarv + 1, :)], prod(size(sample.car)), 1);
            
            tempgfeat(:, end + 1) = reshape(sample.gfeat(1:params.ngfeat, :), prod(size(sample.gfeat)), 1);
            tempbeta(:,:,:,size(tempcam, 2)) = sample.beta;
            
            probrecord2 = [probrecord2, probrecord1(end)];
            temptempsamples = [temptempsamples, sample];
        end
    end
    
    % Disegno la posizione della Camera
    if isempty(tempcam)==0 
        plot3(tempcam(7,:),tempcam(8,:),tempcam(2,:),'o','MarkerFaceColor',Colore,'MarkerEdgeColor',Colore,'MarkerSize',5);    %Get Priors 
    end
        % Disegno la posizione dei TARGET
    if isempty(tempcar)==0
        plot3(tempcar(1,:),tempcar(2,:),tempcar(3,:),'d','MarkerFaceColor',Colore,'MarkerEdgeColor',Colore,'MarkerSize',5);    %Get Priors 
    end
    if isempty(tempper)==0
        plot3(tempper(1,:),tempper(2,:),tempper(3,:),'s','MarkerFaceColor',Colore,'MarkerEdgeColor',Colore,'MarkerSize',5);    %Get Priors 
    end
    
    
    %% N iterations executed for trial
    figure(2);
    subplot(211); plot(probrecord1(500:end));title('Acceptance Ratio 500')
    subplot(212); plot(probrecord2);title('Acceptance Ratio All')
    drawnow
    print('-dpng', ['testfigsmore/samplesfigure' num2str(fn) '_' num2str(trial) '.png']);
    
    if fn > 1
        figure(1);
        print('-dpng', ['testfigsmore/imagefigure' num2str(fn-1), '.png']);
        figure(2);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute MONTECARLO APPROXIMATION for the current trial
    %
    % tempcam,tempcar,tempper...contengono un insieme di Campioni di cui
    % voglio approssimare l'Integrale
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %CAMERA data
    %
    %MEAN VALUE
    Z.cam(:, trial) = mean(tempcam, 2);
    %CAMERA COVARIANCE in Z.V{1}
    Z.V{1, trial} = cov(tempcam');
    %All generated Camera Samples
    Z.cams{trial} = tempcam;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PERSON DATA
    nper = size(sample.per, 2);
    %PERSON MEAN VALUE
    Z.per(:, trial) = mean(tempper, 2);
    %CAVARIANCE in Z.V{2:Number of Person,trial}
    for i = 1:size(sample.per, 2) %get the i-Person       
        try
            disp( [i  (1 + (i-1) * (params.nperv + 1)) i * (params.nperv + 1) - 1] );
            Z.V{i + 1, trial} = cov(tempper((1 + (i-1) * (params.nperv + 1)):(i * (params.nperv + 1) - 1), :)');
            if sample.tcnt(i) == 0
                Z.V{i + 1, trial} = Z.V{i + 1, trial} + diag([0, 0, .3^2, .3^2, 0]);
            end
        catch
            keyboard
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CAR DATA
    ncar = size(sample.car, 2);
    %CAR MEAN VALUE
    Z.car(:, trial) = mean(tempcar, 2);
    %CAVARIANCE in Z.cV
    for i = 1:size(sample.car, 2) % Get the i-Car        
        try
            disp([i (1 + (i-1) * (params.ncarv + 1)) i * (params.ncarv + 1) - 1]);

            Z.cV{i, trial} = cov(tempcar((1 + (i-1) * (params.ncarv + 1)):(i * (params.ncarv + 1) - 1), :)');
            if sample.tccnt(i) == 0
                Z.cV{i, trial} = Z.cV{i, trial} + diag([0, 0, 20^2, 20^2, 0, 0]);
            end
        catch
            keyboard
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% GROUND FEATURES
    Z.gfeat(:, trial) = mean(tempgfeat, 2);
    Z.gfeat(3:3:end, trial) = min(Z.gfeat(3:3:end, trial), 0.99); % 
    Z.gfeat(3:3:end, trial) = max(Z.gfeat(3:3:end, trial), 0.01); % 
    
    for i = 1:size(sample.gfeat, 2)        
        Z.gfV{i, trial} = cov(tempgfeat((1 + (i-1) * params.ngfeat):(i * params.ngfeat - 1), :)') + 0.01^2*eye(2);
    end
    Z.nFeats = size(sample.gfeat, 2);
    
    Z.beta(:,:,:,trial) = mean(tempbeta, 4);
    Z.W(trial) = 1/params.nretry;
%     dbglist
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     SEARCH FOR MAX A POSTERIOR point
    %%%% test max point!
    tsample.cam = Z.cam(:, trial);
    tsample.car = reshape(Z.car(:, trial), 7, length(Z.car(:, trial)) / 7);
    tsample.per = reshape(Z.per(:, trial), 6, length(Z.per(:, trial)) / 6);
    tsample.gfeat = reshape(Z.gfeat(:, trial), 3, length(Z.gfeat(:, trial)) / 3);
    %%POSTERIOR
    [tttprob] =... 
    computeLogPosterior(tsample, prevZ, X, Y, Xc, KLT, corres, corres_car, params);
    display(['[Debug!!!] prob ratio = ' num2str(exp(tttprob - max_prob))]);
    
    Z.W(trial) = (tttprob);
    if (tttprob - max_prob) < log(0.4)
        %UPADATE THE MAX
        Z.cam(:, trial) = maxsample.cam;
        Z.car(:, trial) = reshape(maxsample.car, prod(size(maxsample.car)), 1);
        Z.per(:, trial) = reshape(maxsample.per, prod(size(maxsample.per)), 1);
        Z.gfeat(:, trial) = reshape(maxsample.gfeat, prod(size(maxsample.gfeat)), 1);
        
        Z.gfeat(3:3:end) = min(Z.gfeat(3:3:end), 0.99); % 
        Z.gfeat(3:3:end) = max(Z.gfeat(3:3:end), 0.01); % 
        
        Z.W(trial) = max_prob;
    end
    
    max_prob
    if max_prob == -Inf
        nantrial(end + 1) = trial;
    end
end 
%% All trial Executed: in output genera i campioni per ogni trial
Z.cam(:, nantrial) = [];
Z.car(:, nantrial) = [];
Z.per(:, nantrial) = [];
Z.gfeat(:, nantrial) = [];
Z.V(:, nantrial) = [];
Z.cV(:, nantrial) = [];
Z.gfV(:, nantrial) = [];

Z.tcnt = sample.tcnt + 1;
Z.tccnt = sample.tccnt + 1;

Z.W = Z.W - max(Z.W);
Z.W = exp(Z.W);
Z.W = Z.W ./ sum(Z.W);

Z.gfidx = prevZ.gfidx;
Z.gfcnt = sample.gfcnt + 1;
Z.nSamples = params.nretry - length(nantrial); % cnt - 1;
display(['Average acceptance ratio for one sampling ' num2str(sum(acceptrecord) / sum(samplesrecord), '%.02f')]);
disp(['samples record : ' num2str(samplesrecord)]);
disp(['avgAR = ' num2str(acceptrecord ./ samplesrecord)]);

if sum(acceptrecord) / sum(samplesrecord) < 0.1
%     keyboard;
end
end

%Generate COVARIANCE MATRIX for the PROPOSAL DISTRIBUTION params:
%
% - params.Pert{1} = CAMERA COV. MATRIX
%
% - params.Pert{[2:prevZ.nTargets][prevZ.nTargets:prevZ.nTargets+length(newdet)]}
%                  old detections             new detections 
%
% - params.cPert{[1:prevZ.nCarTargets][prevZ.nCarTargets:length(newcdet)]}
%                  old detections         new detections
%
% - params.PertGF{1:prevZ.nFeats} = GR. FEAT COV
function [params] = getProposalCovariance(prevZ, initidx, newdet, newcdet, params)
%CAMERA COV MATRIX
params.Pert{1} = (prevZ.V{1, initidx} + params.Qcam) ./ params.mPertCam;

%PERSON COV MATRIX
% old targets
for i = 1:prevZ.nTargets
    params.Pert{i + 1} = (prevZ.V{i+1, initidx} + params.Qper2) ./ params.mPertPer;
end
%new detected targets
for i = 1:length(newdet)
    params.Pert{prevZ.nTargets + i + 1} = params.perPert0 / params.mPertPer;
end

%CAR COV MATRIX 
% old targets
for i = 1:prevZ.nCarTargets
    params.cPert{i} = (prevZ.cV{i, initidx} + params.Qcar2) ./ params.mPertCar;
end
%new detected targets
for i = 1:length(newcdet)
    params.cPert{prevZ.nCarTargets + i} = params.carPert0 / params.mPertCar;
end

%GROUND FEATURES COV
for i = 1:prevZ.nFeats
    params.PertGF{i} = (prevZ.gfV{i, initidx}) ./ params.mPertGF;
end

end
% temporary
%%% sample.car(:, j), sample.tccnt(j), tempcar, tempprec, prevZ.W, tempnorm, params, prevZ.nSamples
function prior = mexComputeCarTargetPrior(sample, Z, tid, params)

prior = ones(1, Z.nSamples);

if sample.car((params.ncarv + 1), tid) == 1
    % put height prior always...
    prior = params.normch * exp((sample.car(params.ncarv, tid) - params.mch)^2 / (-2 * params.sch^2)) * prior;
%     prior = params.normch * exp((sample.car(params.ncarv, tid) - params.mch)^2 / (-2 * 0.05^2)) * prior;
else
    prior = params.noCprob * prior;
end

if sample.tccnt(tid) >= 1
    idx = (1 + (tid-1) * (params.ncarv + 1)):(tid * (params.ncarv + 1) - 1);
    tempcar = repmat(sample.car(1:params.ncarv, tid), 1, Z.nSamples) - Z.car(idx, :);
    
    for i = 1:Z.nSamples
        prior(i) = prior(i) * Z.W(i) * (Z.cnormFull(tid, i) * exp(-.5 * diag(tempcar(:, i)' * Z.cprec{tid, i} * tempcar(:, i)))');
    end
    
    % object?
    if sample.car(end, tid) == 1
        prior = prior .* (1e-3 + Z.car(idx(end)+1, :));
    else
        prior = prior .* (1e-3 + (1-Z.car(idx(end)+1, :)));
    end
else
    % nothing
end

end

function [lprob] = computeLogPosterior(sample, prevZ, X, Y, Xc, KLT, corres, corres_car, params)
sample2.per = [params.Aper * sample.per(1:params.nperv, :); sample.per(end, :)];
sample2.car = [params.Acar * sample.car(1:params.ncarv, :); sample.car(end, :)];
sample2.cam = getNextCamera(sample.cam, params);


%% COMPUTE  LIKELYHOOD for OBSERVATION
% no speed defined yet
% lp = 0 - (X.speed - sample2.cam(6))^2/(2*(X.speed/20)^2);
lp = 0;

for i = 1:size(sample.per, 2)
    % compute observation lkhood
    xid = corres(i);
    Xobs = [];
    if xid > 0, Xobs = X.obs(:, xid);  end
    Yobs = [0;0;0]; %Y(:, i);
    lp = lp + mexObservationLkhood(sample2.cam, sample2.per(:, i), Xobs, Yobs, diag(params.Prec1), diag(params.Prec2), params.hdet, params.nhdet, params.lnormDet, params.lnormKtrack);
end

for i = 1:size(sample.car, 2)
    % compute observation lkhood
    xid = corres_car(i);
    Xobs = [];
    if xid > 0, Xobs = Xc.obs(:, xid);  end
    lp = lp + mexCarObservationLkhood(sample2.cam, sample2.car(:, i), Xobs, diag(params.CPrec), params.lnormCDet);
end

for i = 1:size(sample.gfeat, 2)
    % compute observation lkhood
    obs=[KLT.x(i), KLT.y(i)];
    lp = lp + mexFeatureObservationLkhood(sample2.cam, sample.gfeat(:, i), obs, diag(params.kltPrec), params.lkltNorm, params.lkltuprob);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE OVERALL PRIOR

%CAMERA PRIOR
prior(1, :) = mexComputeCompleteCameraPrior(sample.cam, prevZ.cam, prevZ.prec(1, :), prevZ.W, prevZ.normFull(1, :), params.dt);

 if isfield(params, 'vpcam')
    temp = normpdf(sample.cam(4), params.mpcam, params.vpcam);
    prior(1, :) = prior(1, :) * temp;
end
% TARGET PRIOR
for tid = 1:size(sample.per, 2)
        if sample.tcnt(tid) > 0
            tempper = prevZ.per(((tid-1)*(params.nperv + 1)+1):(tid*(params.nperv + 1)), :);
            tempprec = prevZ.prec(tid + 1, :);
            tempnorm = prevZ.normFull(tid+1, :);
        else
            tempper = [];
            tempprec = {};
            tempnorm = [];
        end
        prior(tid + 1, :) = mexComputeTargetPrior(sample.per(:, tid), sample.tcnt(tid), tempper, tempprec, prevZ.W, tempnorm, params, prevZ.nSamples);
end
% CAR PRIOR
for tid = 1:size(sample.car, 2)
    prior_car(tid, :) = mexComputeCarTargetPrior(sample, prevZ, tid, params);
end
%GROUND PRIOR
for tid = 1:size(sample.gfeat, 2)
    prior(size(sample.per, 2) + tid + 1, :) = mexComputeFeturePrior(sample.gfeat(:, tid), sample.gfcnt(tid), prevZ.gfeat(((tid-1)*(params.ngfeat)+1):(tid*(params.ngfeat)), :),...
            prevZ.precgf(tid, :), prevZ.W, prevZ.normgf(tid, :), params, prevZ.nSamples);
end

lprob = lp + log(sum(prod([prior; prior_car], 1)));

end

function [prob] = computePosterior(sample, prevZ, X, Y, Xc, KLT, corres, corres_car, params)
sample2.per = [params.Aper * sample.per(1:params.nperv, :); sample.per(end, :)];
sample2.car = [params.Acar * sample.car(1:params.ncarv, :); sample.car(end, :)];
sample2.cam = getNextCamera(sample.cam, params);

% no speed defined yet
% lp = 0 - (X.speed - sample2.cam(6))^2/(2*(X.speed/20)^2);

for i = 1:size(sample.per, 2)
    % compute observation lkhood
    xid = corres(i);
    Xobs = [];
    if xid > 0, Xobs = X.obs(:, xid);  end
    Yobs = [0;0;0]; %Y(:, i);
    lp = lp + mexObservationLkhood(sample2.cam, sample2.per(:, i), Xobs, Yobs, diag(params.Prec1), diag(params.Prec2), params.hdet, params.nhdet, params.lnormDet, params.lnormKtrack);
end

for i = 1:size(sample.car, 2)
    % compute observation lkhood
    xid = corres_car(i);
    Xobs = [];
    if xid > 0, Xobs = Xc.obs(:, xid);  end
    lp = lp + mexCarObservationLkhood(sample2.cam, sample2.car(:, i), Xobs, diag(params.CPrec), params.lnormCDet);
end

for i = 1:size(sample.gfeat, 2)
    % compute observation lkhood
    obs=[KLT.x(i), KLT.y(i)];
    lp = lp + mexFeatureObservationLkhood(sample2.cam, sample.gfeat(:, i), obs, diag(params.kltPrec), params.lkltNorm, params.lkltuprob);
end

prior(1, :) = mexComputeCompleteCameraPrior(sample.cam, prevZ.cam, prevZ.prec(1, :), prevZ.W, prevZ.normFull(1, :), params.dt);
%%%%%%%%%%%%%%%%%%%%%%%%%% overall prior!
if isfield(params, 'vpcam')
    temp = normpdf(sample.cam(4), params.mpcam, params.vpcam);
    prior(1, :) = prior(1, :) * temp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for tid = 1:size(sample.per, 2)
        if sample.tcnt(tid) > 0
            tempper = prevZ.per(((tid-1)*(params.nperv + 1)+1):(tid*(params.nperv + 1)), :);
            tempprec = prevZ.prec(tid + 1, :);
            tempnorm = prevZ.normFull(tid+1, :);
        else
            tempper = [];
            tempprec = {};
            tempnorm = [];
        end
        prior(tid + 1, :) = mexComputeTargetPrior(sample.per(:, tid), sample.tcnt(tid), tempper, tempprec, prevZ.W, tempnorm, params, prevZ.nSamples);
end

for tid = 1:size(sample.car, 2)
    prior_car(tid, :) = mexComputeCarTargetPrior(sample, prevZ, tid, params);
end

for tid = 1:size(sample.gfeat, 2)
    prior(size(sample.per, 2) + tid + 1, :) = mexComputeFeturePrior(sample.gfeat(:, tid), sample.gfcnt(tid), prevZ.gfeat(((tid-1)*(params.ngfeat)+1):(tid*(params.ngfeat)), :),...
            prevZ.precgf(tid, :), prevZ.W, prevZ.normgf(tid, :), params, prevZ.nSamples);
end

prob = exp(lp) * sum(prod([prior; prior_car], 1));

end

function [prior, prior_car] = initPriors(sample, prevZ, newdet, newcdet, params)

%%%% caching the prior probability for all targets and camera
%CAMERA and PERSONs
prior = zeros(prevZ.nTargets + length(newdet) + 1, prevZ.nSamples);
%CAR
prior_car = zeros(prevZ.nCarTargets + length(newcdet), prevZ.nSamples);

%CAMERA PRIOR
prior(1, :) = mexComputeCompleteCameraPrior(sample.cam, prevZ.cam, prevZ.prec(1, :), prevZ.W, prevZ.normFull(1, :), params.dt);
%%%%%%%%%%%%%%%%%%%%%%%%%% overall prior!
if isfield(params, 'vpcam')
    temp = normpdf(sample.cam(4), params.mpcam, params.vpcam);
    prior(1, :) = prior(1, :) * temp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TARGET PERSON PRIOR
for j = 1:size(sample.per, 2)
    if sample.tcnt(j) > 0
        tempper = prevZ.per(((j-1)*(params.nperv + 1)+1):(j*(params.nperv + 1)), :);
        tempprec = prevZ.prec(j + 1, :);
        tempnorm = prevZ.normFull(j+1, :);
    else
        tempper = [];
        tempprec = {};
        tempnorm = [];
    end
    prior(j + 1, :) = mexComputeTargetPrior(sample.per(:, j), sample.tcnt(j), tempper, tempprec, prevZ.W, tempnorm, params, prevZ.nSamples);
end
%CAR
for j = 1:size(sample.car, 2)
%     if sample.tccnt(j) > 0
%         tempcar = prevZ.car(((j-1)*(params.ncarv + 1)+1):(j*(params.ncarv + 1)), :);
%         tempprec = prevZ.cprec(j, :);
%         tempnorm = prevZ.cnormFull(j, :);
%     else
%         tempcar = [];
%         tempprec = {};
%         tempnorm = [];
%     end
    prior_car(j, :) = mexComputeCarTargetPrior(sample, prevZ, j, params);
end

for j = 1:size(sample.gfeat, 2)
    if sample.gfcnt(j) > 0
        prior(size(sample.per, 2) + j + 1, :) = mexComputeFeturePrior(sample.gfeat(:, j), sample.gfcnt(j), prevZ.gfeat(((j-1)*(params.ngfeat)+1):(j*(params.ngfeat)), :), ...
                prevZ.precgf(j, :), prevZ.W, prevZ.normgf(j, :), params, prevZ.nSamples);
    else
        keyboard
    end
end
end

%function [nsamples, N, Y, Z, prevZ, newdet, newcdet, newgfeats, corres, corres_car] =...
%   initVariables(prevZ, X, Y, Xc, KLT, corres, corres_car, params)
%
% Init the Status of the Variables in the Model Covariance MATRIX of
% CAR/PERSON/CAMERA
% PARAMETERS INPUT:
%
% - prevZ: Estimated Status of Variables in the Model from the previous
%          iteration --> PREDICTION for the current STATE
% - X: PERSON observation      | 
% - Y: Recovered TRACKS        | -->      Current Observation
% - Xc: CAR observation        |
% - KLT: GROUND FEATURES       |
% - corres:  corrispondences for PERSON 
% - corres_car: Corrispondences for Car
% - params:  Model Parameters
%
% PARAMETERS OUTPUT:
%
% - nsamples:  Number of generated samples
% - N:         Total Number of itaration
% - Y:         
% - Z:          Current Status VARIABLE with all Variable unset  
%              
% - prevZ:      OLD Status with updated COVARIANCE MATRIX  (add
%                PERTURBATION)
% - newdet:     new  detections for PEOPLE
% - newcdet:    new  detections for CAR
% - newgfeats:  new  detections for GROUND POINTS
% - corres:     PEOPLE corrispondences updated
% - corres_car: CAR corrispondences updated
function [nsamples, N, Y, Z, prevZ, newdet, newcdet, newgfeats, corres, corres_car] = ...
    initVariables(prevZ, X, Y, Xc, KLT, corres, corres_car, params)
%Number of Samples
nsamples = params.nsamples;

% Total Number of itaration
N = params.burnin + nsamples * params.thinning;

%Generate the CURRENT STATUS VARIABLE with all Variable unset
Z = prevZ;
Z.cam = []; 
Z.per = []; Z.car = [];
Z.gfeat = []; % {x, z, nu} x loc, z loc, valid??
z.gfidx = [];
Z.beta = [];
Z.cams = {};
Z.gfV = {}; %GROUND FEATURE COVARIANCE

% get precision from covariance matrix

%Analyze all the Samples generated in the previous Step:
%For each Sample
%
% 1) Get Covariance Matrix and Precision for GAUSSIAN DISTRIBUTION
for i = 1:prevZ.nSamples
    
    %CAMERA PROCESS: Gaussian Process
    %Get the Covariance Matrix  for the 
    tempMat = prevZ.V{1, i}  + params.Qcam;
    %Get Precision  <--> inverse of covariance MATRIX
    prevZ.prec{1, i} = inv(tempMat);
    %Gaussian Denominator
    prevZ.normFull(1, i) = 1/ sqrt((2 * pi)  ^ size(tempMat, 1) * det(tempMat));
    
    %PERSON PROCESS: Gaussian Process
    %PERSON Variable Distribution
    for j = 1:prevZ.nTargets
        %Covariance 
        tempMat = prevZ.V{j+1, i} + params.Qper2;
        %Precision
        prevZ.prec{j+1, i} = inv(tempMat);
        %Denum
        prevZ.normFull(j + 1, i) = 1/ sqrt((2 * pi)  ^ params.nperv * det(tempMat));
    end
    
    %CAR PROCESS: Gaussian Process
    %CAR Variable Distribution
    for j = 1:prevZ.nCarTargets
        %Covariance
        tempMat = prevZ.cV{j, i} + params.Qcar2;
        %Precision
        prevZ.cprec{j, i} = inv(tempMat);
        %Denum
        prevZ.cnormFull(j, i) = 1/ sqrt((2 * pi)  ^ params.ncarv * det(tempMat));
    end
    
    %GROUND Variable Distribution
    for j = 1:prevZ.nFeats
        prevZ.precgf{j, i} = inv(prevZ.gfV{j, i}); % add some safe guard...to avoid overfitting..
        prevZ.normgf(j, i) = 1/ sqrt((2 * pi)  ^ (params.ngfeat - 1) * det(prevZ.gfV{j, i}));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STORE NEW DETECTIONS
%
%Find if there are new PERSON (DETECTION) in the current frame
newdet = [];
for i = X.idx
    temp = find(i == corres);
    if isempty(temp)
        newdet = [newdet, i]; % Insert the new DETECTION of PERSON
    end
end
%UPDATE the corrispondences for PERSON
corres = [corres, newdet];
%UPDATE the  TARGETS to TRACK with the new DETECTED PERSONS
Y = [Y, zeros(3, length(newdet))];

%GROUND
newgfeats = [];
for i = KLT.idx
    if sum(prevZ.gfidx == i) == 0
        newgfeats(end + 1) = i; %insert new ground features
    end
end

%Find if there are new CAR DETECTION in the current frame
newcdet =[];
for i = Xc.idx
    temp = find(i == corres_car);
    if isempty(temp)
        newcdet = [newcdet, i]; % Insert the new DETECTION of CAR
    end
end
%UPDATE the corrispondences for CAR
corres_car = [corres_car, newcdet];

end
%function [sample] = ...
%getInitialSample(prevZ, Z, X, Xc, KLT, newdet, newcdet, newgfeats, params, initidx)
%
% Start Sampling for the current TIME STEP where:
% - prevZ: is the posterior estimate from the PREVIOUS state
% - Z: is the current STATE VARIABLE
%
% GENERATE SAMPLES FROM THE MOTION MODELS (Linear Gaussian Models)
%    - P(person(t) | person(t-1))
%    - P(car(t ) | car(t-1))
%    - P(ground(t) | ground(t-1))
%
% PARAMETERS INPUT
% - X: PERSON observation
% - Xc: CAR observation
% - Y: TRACKS from previous frame
% - KLT: GROUND FEATURES
% - newdet:  New PERSON DETECTION
% - newcdet: New CAR DETECTION
% - newgfeats: New GROUND DETECTION
% - corres:  Corrispondences for Person
% - corres_car: Corrispondences for Car
% - params:  Model Parameters
% - initidx :   Index of Which Configuration to use
%
% OUTPUT
%
% sample: Sample generated for the VARIABLES in the Models
function [sample] = ...
    getInitialSample(prevZ, Z, X, Xc, KLT, newdet, newcdet, newgfeats, params, initidx)

if isfield(params, 'mpcam');
    mcam = mean(prevZ.cam, 2);
    mcam(4) = params.mpcam;
    sample.cam = mcam;
else
    sample.cam = prevZ.cam(:, initidx);
end
% Get Mean Value from previous state for PERSON and CAR
sample.per = reshape(prevZ.per(:, initidx), (params.nperv + 1), prevZ.nTargets);
sample.car = reshape(prevZ.car(:, initidx), (params.ncarv + 1), prevZ.nCarTargets);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Mean Value  from previous state for GROUND feature states
sample.gfeat = reshape(prevZ.gfeat(:, initidx), (params.ngfeat), prevZ.nFeats);
sample.gfeat(3, :) = rand(1, prevZ.nFeats) < sample.gfeat(3, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERSON - CAR are Multivariate GAUSSIAN PROCESS with:
%
%  mean: Value from previuos state
%  Variance: Variance stored in Z + Additive Noise

%PERSONs
Mper=zeros(params.nperv, 1);
VAddNoise= params.Qper2 + 1e-5 * eye(params.nperv)';
for i = 1:size(sample.per, 2)
    % Old=sample.per(:, i);
    %                                               media nulla                                                  
    sample.per(:, i) = sample.per(:, i) + [mvnrnd(Mper, Z.V{i + 1,initidx} +VAddNoise) 0].';
    
    
    %Disegna
%     figure(11)
%     plot3(Old(1),Old(2),Old(3 ),'og','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',5);
%     plot3(sample.per(1, i),sample.per(2, i),sample.per(3, i),'sg','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',5);
%     linea=[sample.per(1:3, i) Old(1:3)].';
%     line(linea(:,1),linea(:,2),linea(:,3));
end

%CARs
Mcar=zeros(params.ncarv, 1);
VAddNoise= params.Qcar2 + 1e-5 * eye(params.ncarv)';
for i = 1:size(sample.car, 2)
    %Old=sample.car(:, i);
    sample.car(:, i) = sample.car(:, i) + [mvnrnd(Mcar, Z.cV{i, initidx} + VAddNoise) 0].';
    %Disegna
    %figure(11)
    %plot3(Old(1),Old(2),Old(3 ),'or','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5);
    %plot3(sample.car(1, i),sample.car(2, i),sample.car(3, i),'sr','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5);
    %linea=[sample.car(1:3, i) Old(1:3)].';
    %line(linea(:,1),linea(:,2),linea(:,3));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLE INTERACTION STATE Beta(i,j,type) between TARGET (i,j)
%
%Get the number of Samples taking in account the NEW DETECTIONS and the
%number of targets from the previous frame
% 2 --> Indicate the TYPE of interaction modeled ATTRACTION - REPULSION
sample.beta = zeros(prevZ.nTargets + length(newdet), prevZ.nTargets + length(newdet), 2);

%Beta as UPPER TRIANGULAR MATRIX
for i = 1:size(sample.per, 2)
   
    for j = (i+1):size(sample.per,2)
        %Generate the cumulative distribution for beta from the previous
        %state
        ratio = cumsum(prevZ.beta(i, j, :, initidx));
        %Set to 0
        sample.beta(i, j, 1:2) = 0;
        
        %Inverse Sampling
        sample.beta(i, j, sum(rand > ratio) + 1) = 1;
    end
    %A Target cannot' interact with himself
    sample.beta(i, i, 1) = 0; % make square for convenienece: ???
end

%GENERATE SAMPLES FOR NEW DETECTIONS and append to sample variable
%
% PERSON NEW OSERVATION
sample.tcnt = prevZ.tcnt;
for i = newdet
    % compute inverse projection of PERSON
    % randomly sample the humanness by looking at the heights...
    x = X.obs(:, i);
    x = [x(1) + x(3) / 2, x(2) + x(4), x(4)]; %punto a terra
    temp = getIProjection(x', sample.cam);
    
    if sum(isnan(temp)) > 0
        x = x + 2 * randn(1, 3);
        temp = getIProjection(x', sample.cam);
        
    end
    sample.per = [sample.per, [temp; 1]]; % rand < X.pobj(i)]];
    sample.tcnt = [sample.tcnt, 0];
     %Draw the Sample for the New Observation
    % figure(11)
    %plot3(temp(1),temp(2),temp(3),'sr','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5);
    
end

sample.tccnt = prevZ.tccnt;
% CAR NEW OSERVATION (2D) --> ricostruisco il punto a TERRA 3D 
for i = newcdet
    % compute inverse projection of CAR
    % randomly sample the humanness by looking at the heights...
    x = Xc.obs(:, i);
    %PUNTO MEDIO sul piano di terra TERRA
    x = [x(1) + x(3) / 2, x(2) + x(4), x(3), x(4)];
    temp = getCarIProjection(x', sample.cam);
    if sum(isnan(temp)) > 0
        keyboard
    end
    sample.car = [sample.car, [temp; 1]]; % rand < X.pobj(i)]];
    sample.tccnt = [sample.tccnt, 0];
       
    %Draw the Sample for the New Observation
    %figure(11)
    %plot3(temp(1),temp(2),temp(3),'sr','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5);

end

sample.gfcnt = prevZ.gfcnt;
for i = newgfeats
    idx = find(KLT.idx == i);
%         sample.gfeat = [sample.gfeat, [getGIProj([KLT.x(idx); KLT.y(idx)], sample.cam); rand < 0.7]];
    sample.gfeat = [sample.gfeat, [getGIProj([KLT.x(idx); KLT.y(idx)], sample.cam); .7]];
    sample.gfcnt = [sample.gfcnt, 0];
end

%%%%%% initialize interaction state
%     N = size(sample.beta, 2);
for i = 1:length(newdet)
    for j = 1:prevZ.nTargets
        kidx = ceil(2 * rand);
        sample.beta(j, i + prevZ.nTargets, kidx) = 1; % sample!
    end
    for j = 1:(i-1)
        kidx = ceil(2 * rand);
        sample.beta(j + prevZ.nTargets, i + prevZ.nTargets, kidx) = 1; % sample!
    end
    sample.beta(i + prevZ.nTargets, i + prevZ.nTargets, 1) = 0;
end
end

function [ar, prior, prior_car] = getAcceptanceRatio(prevZ, X, Y, Xc, corres, corres_car, KLT, sample, tsample, oid, tid, params, prior, prior_car)
if nargin >= 9
    cachePrior = 1;
else
    error;
    cachePrior = 0;
    prior = [];
    prior_car = [];
end
% compute P(Zt|Xt-1)
% copmute P(Xt|Zmt)
lp1 = 0; % log lkhood/prob for prev sample
lp2 = 0; % log lkhood/prob for current sample

if oid == 0
    % motion model for targets and camera
    sample2.per = [params.Aper * sample.per(1:params.nperv, :); sample.per(end, :)];
    sample2.car = [params.Acar * sample.car(1:params.ncarv, :); sample.car(end, :)];
    sample2.cam = getNextCamera(sample.cam, params);
    
    tsample2.per = [params.Aper * tsample.per(1:params.nperv, :); tsample.per(end, :)];
    tsample2.car = [params.Acar * tsample.car(1:params.ncarv, :); tsample.car(end, :)];
    tsample2.cam = getNextCamera(tsample.cam, params);
    
% no speed defined yet!!!
%     lp1 = lp1 - (X.speed - sample2.cam(6))^2/(2*(X.speed/20)^2);
%     lp2 = lp2 - (X.speed - tsample2.cam(6))^2/(2*(X.speed/20)^2);
    
    for i = 1:size(sample.per, 2)
        % compute observation lkhood
        xid = corres(i);
        Xobs = [];
        if xid > 0, Xobs = X.obs(:, xid);  end

        Yobs = [0;0;0]; %Y(:, i);

        lp1 = lp1 + mexObservationLkhood(sample2.cam, sample2.per(:, i), Xobs, Yobs, diag(params.Prec1), diag(params.Prec2), params.hdet, params.nhdet, params.lnormDet, params.lnormKtrack);
        lp2 = lp2 + mexObservationLkhood(tsample2.cam, tsample2.per(:, i), Xobs, Yobs, diag(params.Prec1), diag(params.Prec2), params.hdet, params.nhdet, params.lnormDet, params.lnormKtrack);
    end
    
    for i = 1:size(sample.car, 2)
        % compute observation lkhood
        xid = corres_car(i);
        Xobs = [];
        if xid > 0, Xobs = Xc.obs(:, xid);  end

        lp1 = lp1 + mexCarObservationLkhood(sample2.cam, sample2.car(:, i), Xobs, diag(params.CPrec), params.lnormCDet);
        lp2 = lp2 + mexCarObservationLkhood(tsample2.cam, tsample2.car(:, i), Xobs, diag(params.CPrec), params.lnormCDet);
    end
    
    for i = 1:size(sample.gfeat, 2)
        % compute observation lkhood
        obs=[KLT.x(i), KLT.y(i)];
        lp1 = lp1 + mexFeatureObservationLkhood(sample2.cam, sample.gfeat(:, i), obs, diag(params.kltPrec), params.lkltNorm, params.lkltuprob) - ...
                mexFeatureObservationLkhood(tsample2.cam, tsample.gfeat(:, i), obs, diag(params.kltPrec), params.lkltNorm, params.lkltuprob);
    end
    
    % no targets' states are changed/won't contribute on the ratio
    lp1 = lp1 - lp2; lp2 = 0;
    if cachePrior 
        lp1 = exp(lp1) * sum(prod([prior; prior_car], 1));
        prior(1, :) = mexComputeCompleteCameraPrior(tsample.cam, prevZ.cam, prevZ.prec(1, :), prevZ.W, prevZ.normFull(1, :), params.dt);
        %%%%%%%%%%%%%%%%%%%%%%%%%% overall prior!
        if isfield(params, 'vpcam')
            temp = normpdf(tsample.cam(4), params.mpcam, params.vpcam);
            prior(1, :) = prior(1, :) * temp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lp2 = sum(prod([prior; prior_car], 1));
    else
        error('');
    end
    % dont need to consider proposal density term since it's gaussian(symmetric)
elseif oid == 1
    % motion model.
    sample2.per = [params.Aper * sample.per(1:params.nperv, tid); sample.per(end, tid)];
    sample2.cam = getNextCamera(sample.cam, params);
    tsample2.per = [params.Aper * tsample.per(1:params.nperv, tid); tsample.per(end, tid)];
    tsample2.cam = getNextCamera(tsample.cam, params);
%     sample2.per = [params.Aper * sample.per(1:params.nperv, :); sample.per(end, :)];
%     tsample2.per = [params.Aper * tsample.per(1:params.nperv, :); tsample.per(end, :)];
    % compute observation lkhood
    xid = corres(tid);
    Xobs = [];
    if xid > 0, Xobs = X.obs(:, xid);  end
    Yobs = [0;0;0]; %Y(:, i);

    lp1 = lp1 + mexObservationLkhood(sample2.cam, sample2.per, Xobs, Yobs, diag(params.Prec1), diag(params.Prec2), params.hdet, params.nhdet, params.lnormDet, params.lnormKtrack);
    lp2 = lp2 + mexObservationLkhood(tsample2.cam, tsample2.per, Xobs, Yobs, diag(params.Prec1), diag(params.Prec2), params.hdet, params.nhdet, params.lnormDet, params.lnormKtrack);
    lp1 = lp1 - lp2; lp2 = 0;
    if cachePrior
%         for i = 1:prevZ.nSamples
%             temp = prevZ.beta(:,:,:,i);
%             idx = find(sample.beta(1:prevZ.nTargets, 1:prevZ.nTargets, :));
%             pbeta1(i) = prod(temp(idx));
%             idx = find(tsample.beta(1:prevZ.nTargets, 1:prevZ.nTargets, :));
%             pbeta2(i) = prod(temp(idx));
%         end
        if params.useInteraction
            [pbeta1, pbeta2] = mexInteractionPrior(sample, tsample, prevZ, params);
        else
            pbeta1 = ones(1, prevZ.nSamples);
            pbeta2 = ones(1, prevZ.nSamples);
        end
%         if (sum(abs(pbeta1 - ttt1) ./ pbeta1) > 1e-10) || (sum(abs(pbeta2 - ttt2) ./ pbeta2) > 1e-10)
%             keyboard;
%         end
            
        lp1 = exp(lp1) * sum(prod([prior; prior_car], 1) .* pbeta1);
%         prior(tid + 1, :) = computeTargetPrior(tsample, prevZ, tid, params);
        if tsample.tcnt(tid) > 0
            tempper = prevZ.per(((tid-1)*(params.nperv + 1)+1):(tid*(params.nperv + 1)), :);
            tempprec = prevZ.prec(tid + 1, :);
            tempnorm = prevZ.normFull(tid+1, :);
        else
            tempper = [];
            tempprec = {};
            tempnorm = [];
        end
        prior(tid + 1, :) = mexComputeTargetPrior(tsample.per(:, tid), tsample.tcnt(tid), tempper, tempprec, prevZ.W, tempnorm, params, prevZ.nSamples);
%         if sum(abs(temp - prior(tid + 1, :)) ./ abs(prior(tid + 1, :))) > 1e-6
%             keyboard
%         end
        lp2 = sum(prod([prior; prior_car], 1) .* pbeta2);
    else
        error('');
    end
    
%     temp = (computePairwisePrior(sample, tid, params)) / (computePairwisePrior(tsample, tid, params));
%     if abs(temp - mexComputeInteraction(sample, tsample, tid, params)) / temp > 1e-10
%         keyboard
%     end
    if params.useInteraction
        lp1 = lp1 * mexComputeInteraction(sample, tsample, tid, params);
    end
elseif oid == 2
    % compute observation lkhood
    sample2.cam = getNextCamera(sample.cam, params);
    tsample2.cam = getNextCamera(tsample.cam, params);
    
    obs=[KLT.x(tid), KLT.y(tid)];
    lp1 = lp1 + mexFeatureObservationLkhood(sample2.cam, sample.gfeat(:, tid), obs, diag(params.kltPrec), params.lkltNorm, params.lkltuprob) - ...
            mexFeatureObservationLkhood(tsample2.cam, tsample.gfeat(:, tid), obs, diag(params.kltPrec), params.lkltNorm, params.lkltuprob);
    if cachePrior
        lp1 = exp(lp1) * sum(prod([prior; prior_car], 1));
        prior(size(sample.per, 2) + tid + 1, :) = mexComputeFeturePrior(tsample.gfeat(:, tid), sample.gfcnt(tid), prevZ.gfeat(((tid-1)*(params.ngfeat)+1):(tid*(params.ngfeat)), :),...
            prevZ.precgf(tid, :), prevZ.W, prevZ.normgf(tid, :), params, prevZ.nSamples);
        lp2 = sum(prod([prior; prior_car], 1));
    else
        error('');
    end
elseif oid == 3
    % motion model.
    sample2.car = [params.Acar * sample.car(1:params.ncarv, tid); sample.car(end, tid)];
    sample2.cam = getNextCamera(sample.cam, params);
    tsample2.car = [params.Acar * tsample.car(1:params.ncarv, tid); tsample.car(end, tid)];
    tsample2.cam = getNextCamera(tsample.cam, params);

    % compute observation lkhood
    xid = corres_car(tid);
    Xobs = [];
    if xid > 0, Xobs = Xc.obs(:, xid);  end

    lp1 = lp1 + mexCarObservationLkhood(sample2.cam, sample2.car, Xobs, diag(params.CPrec), params.lnormCDet);
    lp2 = lp2 + mexCarObservationLkhood(tsample2.cam, tsample2.car, Xobs, diag(params.CPrec), params.lnormCDet);
    lp1 = lp1 - lp2; lp2 = 0;
    
    if cachePrior
        lp1 = exp(lp1) * sum(prod([prior; prior_car], 1));
        prior_car(tid, :) = mexComputeCarTargetPrior(tsample, prevZ, tid, params);
        lp2 = sum(prod([prior; prior_car], 1));
    else
        error('');
    end
end

ar = lp2 / (lp1 + 1e-200);
% if abs(lp2/lp1 - ar) > 0.000001
%     keyboard
% end

end

function Fc = getCameraHesian(cam, params)
Fc = eye(params.ncamv);
if params.ncamv == 8
    Fc(7:8, 5:6) = Fc(7:8, 5:6) + ... 
        [cam(6) * cos(cam(5)) * params.dt, sin(cam(5)) * params.dt; ... 
        -cam(6) * sin(cam(5)) * params.dt, cos(cam(5)) * params.dt];
end
end
% 
% function cam = getNextCamera(cam, params)
% if params.ncamv == 8
%     cam(7) = cam(7) + cam(6) * sin(cam(5)) * params.dt; % + ? - ?
%     cam(8) = cam(8) + cam(6) * cos(cam(5)) * params.dt;
% end
% end

function prob = computePairwisePrior(sample, tid, params)
prob = 1;

for j = 1:size(sample.per, 2)
    if j == tid, continue;  end
    
    % are both an objects??
    if sample.per(6, tid) == 0 || sample.per(6, j) == 0
        prob = prob * params.c1_pair;
    else
        % upper traiangular matrix
        if tid < j
            temp = sample.beta(tid, j, :);
        else
            temp = sample.beta(j, tid, :);
        end
        
        if sum(temp) ~= 1 || sum(temp == 1) == 0
            keyboard
        end
        
        if  temp(1) == 1
            % no interaction
            prob = prob * params.c2_pair;
        elseif temp(2) == 1
            % repulsive interaction
            r = sqrt(sum((sample.per(1:2, tid) - sample.per(1:2, j)).^2));
            prob = prob * ((1/(1+exp(params.i1_pair * (r - params.i2_pair)))) * params.lr1_pair * exp(-1/(params.lr2_pair * r)) + 1e-200);
        elseif temp(3) == 1
            % group movement
            r = (sum((sample.per(3:4, tid) - sample.per(3:4, j)).^2));
            prob = prob * ((1/(1+exp(params.i1_pair * (r - params.i2_pair)))) * params.lg1_pair * exp(-params.lg2_pair * r) + 1e-200);
        else
            error('');
        end
    end
end 
end

function posterior = computeCompletCameraPosterior(sample, prevZ, KLT, params)
% sum_r {P(KLT(t)|Theta(t), Theta(t-1, r), KLT(t-1)) P(Theta(t)|Theta(t-1, r))}
if isempty(KLT.x)
    xx = [];
    yy = [];
    
    posterior = zeros(1, prevZ.nSamples);
    for i = 1:prevZ.nSamples
        camsamples = prevZ.cams{i};
        for j = 1:size(camsamples, 2)
            ncam = getNextCamera(camsamples(:, j), params);
            temp = mvnpdf(sample.cam, ncam, params.Qcam);

            posterior(i) = posterior(i) + temp;
        end
    %     ncam = getNextCamera(prevZ.cam(:, i), params);
    %     prevZ.V{1, i} + params.Qcam
    end
    posterior = posterior ./ size(camsamples, 2);
else
    idx = (KLT.y(:, 1) > min(prevZ.cam(4, :)));
    % filter out non-ground ofs
    xx = KLT.x(idx, :);
    yy = KLT.y(idx, :);
    
    posterior = zeros(1, prevZ.nSamples);
    for i = 1:prevZ.nSamples
        camsamples = prevZ.cams{i};
        for j = 1:size(camsamples, 2)
            tempd = [xx(:, 2),yy(:, 2)]' - getKLTPredictions([xx(:, 1), yy(:, 1)]', camsamples(:, j), sample.cam);
            lkhood = -.5 * diag(tempd' * params.kltPrec * tempd)';
            lkhood(lkhood < -params.klttrunc/2) = -params.klttrunc/2;
            lkhood = sum(lkhood, 2)' * 5 / size(xx, 1);

            ncam = getNextCamera(camsamples(:, j), params);
            temp = mvnpdf(sample.cam, ncam, params.Qcam);

            posterior(i) = posterior(i) + temp;
        end
    %     ncam = getNextCamera(prevZ.cam(:, i), params);
    %     prevZ.V{1, i} + params.Qcam
    end
    posterior = posterior ./ size(camsamples, 2);
end

% 
% tempcam = repmat(sample.cam, 1, prevZ.nSamples) - prevZ.cam;
% for i = 1:prevZ.nSamples
%     ncam = getNextCamera(camsamples(:, j), params);
%     posterior(i) = posterior(i) * prevZ.W(i) * (prevZ.normFull(1, i) * exp(-.5 * diag(tempcam(:, i)' * prevZ.prec{1, i} * tempcam(:, i))));
% end
% % posterior = prior;
% % return;
% if isempty(KLT.x)
%     posterior = prior;
% else
%     idx = (KLT.y(:, 1) > min(prevZ.cam(4, :)));
%     % filter out non-ground ofs
%     xx = KLT.x(idx, :);
%     yy = KLT.y(idx, :);
% 
%     lkhood = zeros(prevZ.nSamples, size(xx, 1));
%     for i = 1:prevZ.nSamples
        cam = prevZ.cam(:, i);
        temp = [xx(:, 2),yy(:, 2)]' - getKLTPredictions([xx(:, 1), yy(:, 1)]', cam, sample.cam);
        lkhood(i, :) = -.5 * diag(temp' * params.kltPrec * temp)'; 
%     end
% 
%     % make the tail of gaussian to be flat -> robust under noise!
%     lkhood(lkhood < -params.klttrunc/2) = -params.klttrunc/2;
% %     lkhood = sum(lkhood, 2)' / size(xx, 1); % do not put too much weight on this???
% %     lkhood = sum(lkhood, 2)' * 2 / size(xx, 1); % do not put too much weight on this???
%     lkhood = sum(lkhood, 2)' * 5 / size(xx, 1); % do not put too much weight on this???
% 
%     % how can we handle numerical problems........
%     posterior = [prior .* (exp(lkhood))];
% end
end

function pred = getKLTPredictions(x, cam1, cam2)

% inverse projection
N = size(x, 2);

z = zeros(2, N);
pred = zeros(2, N);

angle = cam1(5) - cam2(5);
z(2, :) = cam1(1) * cam1(2) ./ (x(2, :) - cam1(4) * ones(1, N));
z(1, :) = (x(1, :) - cam1(3)) .* z(2, :) / cam1(1);
R = [cos(angle), -sin(angle); sin(angle), cos(angle)];
z = R * z;
% projection
pred(1, :) = cam2(1) ./ z(2, :) .* z(1, :) + cam2(3);
pred(2, :) = cam2(1) ./ z(2, :) .* cam2(2) + cam2(4);

end

function prior = computeTargetPrior(sample, Z, tid, params)
prior = ones(1, Z.nSamples);
if sample.per((params.nperv + 1), tid) == 1
    % put height prior always...
    prior = params.normh * exp((sample.per(params.nperv, tid) - params.mh)^2 / (-2 * params.sh^2)) * prior;
    % put velocity prior always...
    prior = params.normv * exp(sum(sample.per(3:4, tid).^2) / (-2 * params.sv^2)) * prior;
else
    prior = params.noHprob * prior;
end

if sample.tcnt(tid) >= 1
    idx = (1 + (tid-1) * (params.nperv + 1)):(tid * (params.nperv + 1) - 1);
    tempper = params.Aper * (repmat(sample.per(1:params.nperv, tid), 1, Z.nSamples) - Z.per(idx, :));
    for i = 1:Z.nSamples
        prior(i) = prior(i) * Z.W(i) * (Z.normFull(tid + 1, i) * exp(-.5 * diag(tempper(:, i)' * Z.prec{tid + 1, i} * tempper(:, i)))');
    end
    % object?
    if sample.per(end, tid) == 1
        prior = prior .* Z.per(idx(end)+1, :);
    else
        prior = prior .* (1-Z.per(idx(end)+1, :));
    end
else
    % nothing
end

end

function lkhood = getObservationLKhood(cam, state, Xobs, Yobs, params)
lkhood = 0;
ImPred = getProjection(state, cam);

if ~isempty(Xobs) % is there observation(detection)?
    ImObs = Xobs;
    ImObs = [ImObs(1) + ImObs(3) / 2; ImObs(2) + ImObs(4); ImObs(4)];
    temp = ImObs - ImPred';
    lkhood = lkhood + -.5 * temp' * (params.Prec1 / ((ImObs(3)/128) ^ 2)) * temp + params.lnormDet;
    
    % P(nodet|human)
    if state(6) == 1
        lkhood = lkhood + log(params.hdet);
    else
        lkhood = lkhood + log(params.nhdet);
    end
else
    % P(nodet|human)
    if state(6) == 1
        lkhood = log(1 - params.hdet);
    else
        lkhood = log(1 - params.nhdet);
    end
end

if sum(Yobs) ~= 0 % is there valid ms track?
    Yobs(2) = Yobs(2) + Yobs(3) / 2;
    temp = Yobs - ImPred';
    lkhood = lkhood + -.5 * temp' * (params.Prec2 / ((Yobs(3)/128)^2)) * temp + params.lnormKtrack;
end

end

function h = computeEntrophy(z, idx, rKernel)
% h = 0;
% 
% ndim = size(z.state, 1);
% iKernel = inv(rKernel);
% 
% for i = idx
%     temph = 0;
%     for j = idx
%         tempd = z.state(:, i) - z.state(:, j);
%         temph = temph + z.w(j) * exp(-.5 * tempd * iKernel * tempd);
%     end
%     
%     h = h + z.w(i) * log(temph /(sqrt(2 * pi) * rKernel)^ndim);
% end
% 
% h = -h;
end