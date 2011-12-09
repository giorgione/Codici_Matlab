function [Z, corres] = MCMCSamplesJointStates(prevZ, X, Y, corres, params)
% compute uncertainty for all targets
% just use parameter value for now..
nsamples = params.nsamples;
N = params.burnin + nsamples * params.thinning;

% compute sampling weight for each target and camera parameters
% tweight = ones(1, prevZ.nTargets + 1) / (prevZ.nTargets + 1);
% tweight = cumsum(tweight);

cnt = 1;
Z = prevZ;
Z.cam = [];
Z.per = [];

% precompute normalization parameters!!!!
params.Prec1 = inv(params.Qdet); % observation noise
params.Prec2 = inv(params.Qktrack); % ms track noise.
params.PrecCam = inv(params.Qcam); 
params.PrecPer1 = inv(params.Aper * params.Qper1 * params.Aper'); % inv(params.Qper);
params.PrecPer2 = inv(params.Aper * params.Qper2 * params.Aper'); % inv(params.Qper);
% params.PrecPer1 = inv(params.Qper1); % inv(params.Qper);
% params.PrecPer2 = inv(params.Qper2); % inv(params.Qper);

params.normh = 1/sqrt(2 * pi * 0.1^2);
params.normFull1 = 1/ sqrt((2 * pi)  ^ params.nperv * det(params.Aper * params.Qper1 * params.Aper'));
params.normFull2 = 1/ sqrt((2 * pi)  ^ params.nperv * det(params.Aper * params.Qper2 * params.Aper'));
% params.normFull1 = 1/ sqrt((2 * pi)  ^ params.nperv * det(params.Qper1));
% params.normFull2 = 1/ sqrt((2 * pi)  ^ params.nperv * det(params.Qper2));
params.lnormDet = log(1/ sqrt((2 * pi)  ^ 3 * det(params.Qdet)));
params.lnormKtrack = log(1/ sqrt((2 * pi)  ^ 3 * det(params.Qktrack)));

params.hdet = 0.6;
params.nhdet = 0.1; 

params.mh = 1.7;
params.sh = 0.1;
params.noHprob = 0.06;

params.flipObj = 0.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
newdet = [];
for i = X.idx
    temp = find(i == corres);
    if isempty(temp)
        newdet = [newdet, i];
    end
end
corres = [corres, newdet];
Y = [Y, zeros(3, length(newdet))];

samplesrecord = zeros(1, prevZ.nTargets + 1 + length(newdet));
% initialize the sample
% sample from prevZ
% maybe we can go for multiple times to guarantee multi-modality
for trial = 1:params.nretry
    initidx = ceil(rand * prevZ.nSamples);

    sample.cam = prevZ.cam(:, initidx);
    sample.per = reshape(prevZ.per(:, initidx), (params.nperv + 1), prevZ.nTargets);
    sample.tcnt = prevZ.tcnt;
    for i = newdet
        % compute inverse projection
        % randomly sample the humanness by looking at the heights...
        x = X.obs(:, i);
        x = [x(1) + x(3) / 2, x(2) + x(4), x(4)];
        sample.per = [sample.per, [getIProjection(x', sample.cam); rand < X.pobj(i)]];
        sample.tcnt = [sample.tcnt, 0];
    end
    tempcnt = 0;
    %%%%%%%%%%%%%%%%%
%     tempdbg = zeros(prevZ.nTargets * params.nperv, N);
%     dbglist = zeros(1, N);
    %%%%%%%%%%%%
    %%%% caching the prior probability for all targets and camera
    prior = zeros(prevZ.nTargets + length(newdet) + 1, prevZ.nSamples);
    tempcam = repmat(sample.cam, 1, prevZ.nSamples) - prevZ.cam;
    prior(1, :) = exp(-.5 * diag(tempcam' * params.PrecCam * tempcam))';
    
    tsample = sample;
%     tsample.per(1:params.nperv, :) = params.Aper * tsample.per(1:params.nperv, :);
    for j = 1:size(sample.per, 2)
        prior(j + 1, :) = computeTargetPrior(tsample, prevZ, j, params);
    end
    
    for i = 1:N
        while(1)
            tempcnt = tempcnt + 1;
            % get sample
            % select a target or camera
            tid = ceil(rand * (prevZ.nTargets + length(newdet) + 1));

            tsample = sample;
            if tid == 1 % camera
                % need to be able to deal with multimodality
                tsample.cam = tsample.cam + mvnrnd(zeros(1, params.ncamv), params.camPert)';
            else
                if sample.tcnt(tid-1) > 1
                    temp = tsample.per(end, tid-1);
                    if rand < params.flipObj
                        temp = ~temp;
                    end
                    tsample.per(:, tid-1) = [tsample.per(1:params.nperv, tid-1) + mvnrnd(zeros(1, params.nperv), params.perPert2)'; temp];
                elseif sample.tcnt(tid-1) == 1
                    temp = tsample.per(end, tid-1);
                    if rand < params.flipObj
                        temp = ~temp;
                    end
                    tsample.per(:, tid-1) = [tsample.per(1:params.nperv, tid-1) + mvnrnd(zeros(1, params.nperv), params.perPert1)'; temp];
                else
                    temp = tsample.per(end, tid-1);
                    if rand < params.flipObj
                        temp = ~temp;
                    end
                    tsample.per(:, tid-1) = [tsample.per(1:params.nperv, tid-1) + mvnrnd(zeros(1, params.nperv), params.perPert0)'; temp];
                end
            end

            % evaluate acceptance ratio
            [ar, tprior] = getAcceptanceRatio(prevZ, X, Y, corres, sample, tsample, tid, params, prior); % compute acceptance ratio
%             if isnan(ar)
%                 keyboard
%             end
%             if abs(ar - getAcceptanceRatio(prevZ, X, Y, corres, sample, tsample, tid, params))/ar > 1e-3
%                 error('fatal!');
%             end
            % discard/accept
%             dbglist(i) = dbglist(i) + 1;
            if ar > rand, break; end
        end
        % accepted sample
        sample = tsample;  
        prior = tprior;
    %     tempcnt
    %     tsample.cam
        %%%%%%%%%%
%         tempdbg(:, i) = reshape(sample.per, prevZ.nTargets * params.nperv, 1);
        %%%%%%%%%%
        if i > params.burnin
            if mod(i - params.burnin, params.thinning) ~= 0, continue; end
            % save them.
            Z.cam(:, cnt) = sample.cam;
            Z.per(:, cnt) = reshape([params.Aper * sample.per(1:params.nperv, :); sample.per(params.nperv + 1, :)], prod(size(sample.per)), 1);
            cnt = cnt + 1;
        end
        
        samplesrecord(tid) = samplesrecord(tid) + 1;
    end
%     dbglist
end

Z.tcnt = sample.tcnt + 1;
Z.nSamples = cnt - 1;
display(['Average # of trial for one sampling ' num2str(tempcnt/N, '%.02f')]);
samplesrecord
end

function [ar, prior] = getAcceptanceRatio(prevZ, X, Y, corres, sample, tsample, tid, params, prior)
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

if tid == 1
    % motion model.
    sample2.per = [params.Aper * sample.per(1:params.nperv, :); sample.per(end, :)];
    tsample2.per = [params.Aper * tsample.per(1:params.nperv, :); tsample.per(end, :)];
    
    for i = 1:size(sample.per, 2)
        % compute observation lkhood
        xid = corres(i);
        Xobs = [];
        if xid > 0, Xobs = X.obs(:, xid);  end
        Yobs = Y(:, i);
        lp1 = lp1 + getObservationLKhood(sample.cam, sample2.per(:, i), Xobs, Yobs, params);
        lp2 = lp2 + getObservationLKhood(tsample.cam, tsample2.per(:, i), Xobs, Yobs, params);
    end
    % no targets' states are changed/won't contribute on the ratio
    lp1 = lp1 - lp2; lp2 = 0;
    if cachePrior 
        lp1 = exp(lp1) * sum(prod(prior, 1));
        tempcam = repmat(tsample.cam, 1, prevZ.nSamples) - prevZ.cam;
        prior(1, :) = exp(-.5 * diag(tempcam' * params.PrecCam * tempcam))';
        lp2 = sum(prod(prior, 1));
    else
    end
    % dont need to consider proposal density term since it's gaussian(symmetric)
else
    tid = tid - 1;
    % motion model.
    sample2.per = [params.Aper * sample.per(1:params.nperv, :); sample.per(end, :)];
    tsample2.per = [params.Aper * tsample.per(1:params.nperv, :); tsample.per(end, :)];
    % compute observation lkhood
    xid = corres(tid);
    Xobs = [];
    if xid > 0, Xobs = X.obs(:, xid);  end
    Yobs = Y(:, tid);
    lp1 = lp1 + getObservationLKhood(sample.cam, sample2.per(:, tid), Xobs, Yobs, params);
    lp2 = lp2 + getObservationLKhood(tsample.cam, tsample2.per(:, tid), Xobs, Yobs, params);
    
    lp1 = lp1 - lp2; lp2 = 0;
    if cachePrior 
        lp1 = exp(lp1) * sum(prod(prior, 1));
        prior(tid + 1, :) = computeTargetPrior(tsample, prevZ, tid, params);
        lp2 = sum(prod(prior, 1));
    else
    end
end

ar = lp2 / lp1;

end

function prior = computeTargetPrior(sample, Z, tid, params)
prior = ones(1, Z.nSamples);
if sample.per((params.nperv + 1), tid) == 1
    prior = params.normh * exp((sample.per(params.nperv, tid) - params.mh)^2 / (-2 * params.sh^2)) * prior;
else
    prior = params.noHprob * prior;
end
if sample.tcnt(tid) > 1
    idx = (1 + (tid-1) * (params.nperv + 1)):(tid * (params.nperv + 1) - 1);
    tempper = params.Aper * (repmat(sample.per(1:params.nperv, tid), 1, Z.nSamples) - Z.per(idx, :));
    prior = prior .* (params.normFull2 * exp(-.5 * diag(tempper' * params.PrecPer2 * tempper))');
    prior = prior .* (0.1 .^ ((sample.per(end, tid) * ones(1, Z.nSamples) ~= Z.per(idx(end)+1, :))));
elseif sample.tcnt(tid) == 1
    idx = (1 + (tid-1) * (params.nperv + 1)):(tid * (params.nperv + 1) - 1);
    tempper = params.Aper * (repmat(sample.per(1:params.nperv, tid), 1, Z.nSamples) - Z.per(idx, :));
    prior = prior .* (params.normFull1 * exp(-.5 * diag(tempper' * params.PrecPer1 * tempper))');
    prior = prior .* (0.1 .^ ((sample.per(end, tid) * ones(1, Z.nSamples) ~= Z.per(idx(end)+1, :))));
else
    % nothing
%     if sample.per((params.nperv + 1), tid) == 1
%         prior = params.normh * exp((sample.per(params.nperv, tid) - params.mh)^2 / (-2 * params.sh^2))' * ones(1, Z.nSamples);
%     else
%         prior = params.noHprob * ones(1, Z.nSamples);
%     end
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