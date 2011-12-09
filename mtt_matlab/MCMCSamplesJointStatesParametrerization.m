function [Z, corres] = MCMCSamplesJointStatesParametrerization(prevZ, X, Xc, Y, KLT, corres, params)
nsamples = params.nsamples;
N = params.burnin + nsamples * params.thinning;

% cnt = 1;
Z = prevZ; Z.cam = []; Z.per = [];
Z.gfeat = []; % {x, z, nu} x loc, z loc, valid??
z.gfidx = []; Z.beta = []; Z.cams = {}; Z.gfV = {};

% get precision from covariance matrix
for i = 1:prevZ.nSamples
    %%%%%%%%% CAMMOVEMENT!!!!!!!!
%     Fc = getCameraHesian(prevZ.cam(:, i), params);
%     tempMat = Fc * (prevZ.V{1, i}  + params.Qcam) * Fc';
    tempMat = prevZ.V{1, i}  + params.Qcam;
%     tempMat = params.Qcam;
    prevZ.prec{1, i} = inv(tempMat);
    
    prevZ.normFull(1, i) = 1/ sqrt((2 * pi)  ^ size(tempMat, 1) * det(tempMat));
    for j = 1:prevZ.nTargets
%         tempMat = params.Aper * (prevZ.V{j+1, i} + params.Qper2) * params.Aper';
        tempMat = prevZ.V{j+1, i} + params.Qper2;
        prevZ.prec{j+1, i} = inv(tempMat);
        prevZ.normFull(j + 1, i) = 1/ sqrt((2 * pi)  ^ params.nperv * det(tempMat));
    end
    
    for j = 1:prevZ.nCarTargets
%         tempMat = params.Aper * (prevZ.V{j+1, i} + params.Qper2) * params.Aper';
        tempMat = prevZ.cV{j, i} + params.Qcar2;
        prevZ.cprec{j, i} = inv(tempMat);
        prevZ.cnormFull(j, i) = 1/ sqrt((2 * pi)  ^ params.ncarv * det(tempMat));
    end
        
    for j = 1:prevZ.nFeats
        prevZ.precgf{j, i} = inv(prevZ.gfV{j, i}); % add some safe guard...to avoid overfitting..
        prevZ.normgf(j, i) = 1/ sqrt((2 * pi)  ^ (params.ngfeat - 1) * det(prevZ.gfV{j, i}));
    end
end

newdet = [];
for i = X.idx
    temp = find(i == corres);
    if isempty(temp)
        newdet = [newdet, i];
    end
end
corres = [corres, newdet];
Y = [Y, zeros(3, length(newdet))];

newgfeats = [];
for i = KLT.idx
    if sum(prevZ.gfidx == i) == 0
        newgfeats(end + 1) = i;
    end
end

samplesrecord = zeros(1, prevZ.nTargets + 1 + length(newdet) + prevZ.nFeats);
acceptrecord = zeros(1, prevZ.nTargets + 1 + length(newdet) + prevZ.nFeats);
% initialize the sample
% sample from prevZ
% maybe we can go for multiple times to guarantee multi-modality
for trial = 1:params.nretry
    initidx = ceil(rand * prevZ.nSamples);
    if isfield(params, 'mpcam');
        mcam = mean(prevZ.cam, 2);
        
        mcam(4) = params.mpcam;
        
        sample.cam = mcam;
%         sample = getMCMCInitialization(X, mcam, params); % prevZ.cam(:, initidx);
%         mcam'
%         sample.cam'
    else
        sample.cam = prevZ.cam(:, initidx);
    end
    sample.per = reshape(prevZ.per(:, initidx), (params.nperv + 1), prevZ.nTargets);
    %%%%% ground feature states
    sample.gfeat = reshape(prevZ.gfeat(:, initidx), (params.ngfeat), prevZ.nFeats);
    sample.gfeat(3, :) = rand(1, prevZ.nFeats) < sample.gfeat(3, :);
    %%%%%% add uncertainty??
%     sample.cam = sample.cam; % + mvnrnd(zeros(params.ncamv, 1), Z.V{1,initidx})';
%     sample.cam = getNextCamera(sample.cam, params); % pretranslate camera!
    
    for i = 1:size(sample.per, 2)
        sample.per(:, i) = sample.per(:, i) + [mvnrnd(zeros(5, 1), Z.V{i + 1,initidx} + params.Qper2 + 1e-5 * eye(5))'; 0];
        % initialization trick!
%         sample.per(5, i) = 1.7 + randn * 0.1;
    end
    %%%%%% sample interaction state
    sample.beta = zeros(prevZ.nTargets + length(newdet), prevZ.nTargets + length(newdet), 2);
    for i = 1:size(sample.per, 2)
        for j = (i+1):size(sample.per,2)
            ratio = cumsum(prevZ.beta(i, j, :, initidx));
            sample.beta(i, j, 1:2) = 0;
            sample.beta(i, j, sum(rand > ratio) + 1) = 1;
        end
        sample.beta(i, i, 1) = 0; % make square for convenienece
    end
        
    sample.tcnt = prevZ.tcnt;
    for i = newdet
        % compute inverse projection
        % randomly sample the humanness by looking at the heights...
        x = X.obs(:, i);
        x = [x(1) + x(3) / 2, x(2) + x(4), x(4)];
        temp = getIProjection(x', sample.cam);
        if sum(isnan(temp)) > 0
            x = x + 2 * randn(1, 3);
            temp = getIProjection(x', sample.cam);
        end
        sample.per = [sample.per, [temp; 1]]; % rand < X.pobj(i)]];
        sample.tcnt = [sample.tcnt, 0];
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
    tempcnt = 0;
    
    %%%% caching the prior probability for all targets and camera
    prior = zeros(prevZ.nTargets + length(newdet) + 1, prevZ.nSamples);
    prior(1, :) = mexComputeCompleteCameraPrior(sample.cam, prevZ.cam, prevZ.prec(1, :), prevZ.W, prevZ.normFull(1, :), params.dt);
    %%%%%%%%%%%%%%%%%%%%%%%%%% overall prior!
    if isfield(params, 'vpcam')
        temp = normpdf(sample.cam(4), params.mpcam, params.vpcam);
        prior(1, :) = prior(1, :) * temp;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    for j = 1:size(sample.gfeat, 2)
        if sample.gfcnt(j) > 0
            prior(size(sample.per, 2) + j + 1, :) = mexComputeFeturePrior(sample.gfeat(:, j), sample.gfcnt(j), prevZ.gfeat(((j-1)*(params.ngfeat)+1):(j*(params.ngfeat)), :), ...
                    prevZ.precgf(j, :), prevZ.W, prevZ.normgf(j, :), params, prevZ.nSamples);
        else
            keyboard
        end
    end
    
    if (sum(sum(isnan(prior))) > 0) || (sum(sum(prior == inf)) > 0) || (sum(sum(prior <= 0)) > 0)
        keyboard;
    end
%     prior
    %%%% get the perturbation matrices relative to the variation..
    params.Pert{1} = (prevZ.V{1, initidx} + params.Qcam) ./ params.mPertCam;
%     params.Pert{1} = (params.Qcam) ./ params.mPertCam;
%     params.Pert{1} = diag(min(diag(params.Pert{1})', [1, 0.01, 1e-10, 1, 0.2*pi/180, 0.05, 1e-10, 1e-10].^2));
    for i = 1:Z.nTargets
        params.Pert{i + 1} = (prevZ.V{i+1, initidx} + params.Qper2) ./ params.mPertPer;
    end
    for i = 1:length(newdet)
        params.Pert{Z.nTargets + i + 1} = params.perPert0 / params.mPertPer;
    end
    for i = 1:Z.nFeats
        params.PertGF{i} = (prevZ.gfV{i, initidx}) ./ params.mPertGF;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tempcam = [];
    tempper = [];
    tempbeta = [];
    tempgfeat = [];
    
    if isfield(params, 'camscnt')
        perTargSample = [params.camscnt * params.thinning, params.perscnt * params.thinning * ones(1, prevZ.nTargets + length(newdet)), params.featscnt * params.thinning * ones(1, prevZ.nFeats)];
        cumSamples = cumsum(perTargSample) / sum(perTargSample);
        N = params.burnin + sum(perTargSample);
    end
    
    for i = 1:N
        dbgcnt = 0;
        tempcnt = tempcnt + 1;
        % get sample
        % select a target or camera
        if isfield(params, 'camscnt')
            tid = ceil(rand * (prevZ.nFeats + prevZ.nTargets + length(newdet) + 1));
            tid = sum(rand > cumSamples) + 1;
        else
            tid = ceil(rand * (prevZ.nFeats + prevZ.nTargets + length(newdet) + 1));
        end
        
        oid = 0; tsample = sample;
        sid = 0;
        if tid == 1
            oid = 0; sid = tid;
            tsample.cam = tsample.cam + mvnrnd(zeros(1, params.ncamv), params.Pert{1})';
            if (params.useCamAnnealing == 1)
                params.Pert{1} = params.Pert{1} .* params.AnnealingConst;
            end
        elseif tid <= 1+size(sample.per, 2)
            oid = 1; sid = tid - 1;
            temp = tsample.per(end, tid-1);
            if rand < params.flipObj
                temp = ~temp;
            end
            tsample.per(:, tid-1) = [tsample.per(1:params.nperv, tid-1) + mvnrnd(zeros(1, params.nperv), params.Pert{tid})'; temp];
            % sample beta!

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
        else
            oid = 2; sid = tid - 1 - size(sample.per, 2);
            temp = tsample.gfeat(end, sid);
            if rand < params.flipFeat
                temp = ~temp;
            end
            tsample.gfeat(:, sid) = [tsample.gfeat(1:(params.ngfeat-1), sid) + mvnrnd(zeros(1, params.ngfeat-1), params.PertGF{sid})'; temp];
        end
        
        % evaluate acceptance ratio
        [ar, tprior] = getAcceptanceRatio(prevZ, X, Y, corres, KLT, sample, tsample, oid, sid, params, prior); % compute acceptance ratio
        if isnan(ar)
%             keyboard
        end

        % accepted sample
        if ar > rand
            sample = tsample;  
            prior = tprior;
            acceptrecord(tid) = acceptrecord(tid) + 1;
        end
        samplesrecord(tid) = samplesrecord(tid) + 1;
        
        if i > params.burnin
            if mod(i - params.burnin, params.thinning) ~= 0, continue; end
            % save them.
            tempcam(:, end + 1) = getNextCamera(sample.cam, params);
            tempper(:, end + 1) = reshape([params.Aper * sample.per(1:params.nperv, :); sample.per(params.nperv + 1, :)], prod(size(sample.per)), 1);
            tempgfeat(:, end + 1) = reshape(sample.gfeat(1:params.ngfeat, :), prod(size(sample.gfeat)), 1);
            tempbeta(:,:,:,size(tempcam, 2)) = sample.beta;
%             Z.cam(:, cnt) = sample.cam;
%             Z.per(:, cnt) = reshape([params.Aper * sample.per(1:params.nperv, :); sample.per(params.nperv + 1, :)], prod(size(sample.per)), 1);
%             cnt = cnt + 1;
        end
    end
    
%     prior
    Z.cam(:, trial) = mean(tempcam, 2);
    Z.V{1, trial} = cov(tempcam');
    Z.cams{trial} = tempcam;
    
    nper = size(sample.per, 2);
    Z.per(:, trial) = mean(tempper, 2);
    for i = 1:size(sample.per, 2)        
        try
            Z.V{i + 1, trial} = cov(tempper((1 + (i-1) * (params.nperv + 1)):(i * (params.nperv + 1) - 1), :)');
            if sample.tcnt(i) == 0
                Z.V{i + 1, trial} = Z.V{i + 1, trial} + diag([0, 0, .3^2, .3^2, 0]);
            end
        catch
            keyboard
        end
    end
    
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
end

Z.tcnt = sample.tcnt + 1;
Z.gfidx = prevZ.gfidx;
Z.gfcnt = sample.gfcnt + 1;
Z.nSamples = params.nretry; % cnt - 1;
display(['Average acceptance ratio for one sampling ' num2str(sum(acceptrecord) / sum(samplesrecord), '%.02f')]);
disp(['samples record : ' num2str(samplesrecord)]);
disp(['avgAR = ' num2str(acceptrecord ./ samplesrecord)]);

if sum(acceptrecord) / sum(samplesrecord) < 0.1
%     keyboard;
end
end

function [ar, prior] = getAcceptanceRatio(prevZ, X, Y, corres, KLT, sample, tsample, oid, tid, params, prior)
if nargin >= 9
    cachePrior = 1;
else
    error;
    cachePrior = 0;
    prior = [];
end
% compute P(Zt|Xt-1)
% copmute P(Xt|Zmt)
lp1 = 0; % log lkhood/prob for prev sample
lp2 = 0; % log lkhood/prob for current sample

if oid == 0
    % motion model for targets and camera
    sample2.per = [params.Aper * sample.per(1:params.nperv, :); sample.per(end, :)];
    sample2.cam = getNextCamera(sample.cam, params);
    tsample2.per = [params.Aper * tsample.per(1:params.nperv, :); tsample.per(end, :)];
    tsample2.cam = getNextCamera(tsample.cam, params);
    
    lp1 = lp1 - (X.speed - sample2.cam(6))^2/(2*(X.speed/20)^2);
    lp2 = lp2 - (X.speed - tsample2.cam(6))^2/(2*(X.speed/20)^2);
    
    for i = 1:size(sample.per, 2)
        % compute observation lkhood
        xid = corres(i);
        Xobs = [];
        if xid > 0, Xobs = X.obs(:, xid);  end
        Yobs = Y(:, i);

        lp1 = lp1 + mexObservationLkhood(sample2.cam, sample2.per(:, i), Xobs, Yobs, diag(params.Prec1), diag(params.Prec2), params.hdet, params.nhdet, params.lnormDet, params.lnormKtrack);
        lp2 = lp2 + mexObservationLkhood(tsample2.cam, tsample2.per(:, i), Xobs, Yobs, diag(params.Prec1), diag(params.Prec2), params.hdet, params.nhdet, params.lnormDet, params.lnormKtrack);
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
        lp1 = exp(lp1) * sum(prod(prior, 1));
        prior(1, :) = mexComputeCompleteCameraPrior(tsample.cam, prevZ.cam, prevZ.prec(1, :), prevZ.W, prevZ.normFull(1, :), params.dt);
        %%%%%%%%%%%%%%%%%%%%%%%%%% overall prior!
        if isfield(params, 'vpcam')
            temp = normpdf(tsample.cam(4), params.mpcam, params.vpcam);
            prior(1, :) = prior(1, :) * temp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lp2 = sum(prod(prior, 1));
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
    Yobs = Y(:, tid);

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
            
        lp1 = exp(lp1) * sum(prod(prior, 1) .* pbeta1);
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
        lp2 = sum(prod(prior, 1) .* pbeta2);
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
        lp1 = exp(lp1) * sum(prod(prior, 1));
        prior(size(sample.per, 2) + tid + 1, :) = mexComputeFeturePrior(tsample.gfeat(:, tid), sample.gfcnt(tid), prevZ.gfeat(((tid-1)*(params.ngfeat)+1):(tid*(params.ngfeat)), :),...
            prevZ.precgf(tid, :), prevZ.W, prevZ.normgf(tid, :), params, prevZ.nSamples);
        lp2 = sum(prod(prior, 1));
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