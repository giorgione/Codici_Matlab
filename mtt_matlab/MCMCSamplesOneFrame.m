function [Z] = MCMCSamplesOneFrame(prevZ, X, params)
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

params.Prec1 = inv(params.Qdet); % observation noise
params.PrecCam = inv(params.Qcam); 

nTarget = size(X.obs, 2);
samplesrecord = zeros(1, nTarget + 1);

for trial = 1:params.nretry
    initidx = ceil(rand * prevZ.nSamples);

    sample.cam = prevZ.cam(:, initidx);
    for i = 1:nTarget
        % compute inverse projection
        % randomly sample the humanness by looking at the heights...
        x = X.obs(:, i);
        x = [x(1) + x(3) / 2, x(2) + x(4), x(4)];
        sample.per(:, i) = [getIProjection(x', sample.cam); rand < X.pobj(i)];
    end

    %%%% caching the prior probability for all targets and camera
    prior = zeros(nTarget + 1, prevZ.nSamples);
    tempcam = repmat(sample.cam, 1, prevZ.nSamples) - prevZ.cam;
    prior(1, :) = exp(-.5 * diag(tempcam' * params.PrecCam * tempcam))';
    for j = 1:size(sample.per, 2)
        if sample.per(6, j) == 1
            prior(j + 1, :) = 1/sqrt(2 * pi * 0.1^2) * exp((sample.per(5, j)-1.7)^2/(-2 * 0.1^2))' * ones(1, prevZ.nSamples); % * X.pobj(i);
        else
            prior(j + 1, :) = 1/3 * ones(1, prevZ.nSamples); % * (1-X.pobj(i));
        end
    end
    
    tempcnt = 0;
    for i = 1:N
        tempcnt2 = 0;
        while(1)
            tempcnt = tempcnt + 1;
            tempcnt2 = tempcnt2 + 1;
            % get sample
            % select a target or camera
            tid = ceil(rand * (nTarget + 1));

            tsample = sample;
            if tid == 1 % camera
                % need to be able to deal with multimodality
                tsample.cam = tsample.cam + mvnrnd(zeros(1, params.ncamv), params.camPert)';
            else
                temp = tsample.per(end, tid-1);
                if rand < 0.2
                    temp = ~temp;
                end
                tsample.per(:, tid-1) = [tsample.per(1:5, tid-1) + mvnrnd(zeros(1, params.nperv), params.perPert)'; temp];
            end

            % evaluate acceptance ratio
            [ar, tprior] = getAcceptanceRatio(prevZ, X, sample, tsample, tid, params, prior); % compute acceptance ratio
%             if isnan(ar)
%                 keyboard
%             end
%             if abs(ar - getAcceptanceRatio(prevZ, X, Y, corres, sample, tsample, tid, params))/ar > 1e-3
%                 error('fatal!');
%             end
            % discard/accept
%             dbglist(i) = dbglist(i) + 1;
            if ar > rand, break; end
            if tempcnt2 > 200, break; end;
%             break;
        end
        if tempcnt2 > 200, break; end;
        % accepted sample
%         if ar > rand
            sample = tsample;  
            prior = tprior;
%         end            
    %     tempcnt
    %     tsample.cam
        %%%%%%%%%%
%         tempdbg(:, i) = reshape(sample.per, prevZ.nTargets * params.nperv, 1);
        %%%%%%%%%%
        if i > params.burnin
            if mod(i - params.burnin, params.thinning) ~= 0, continue; end
            % save them.
            Z.cam(:, cnt) = sample.cam;
            Z.per(:, cnt) = reshape(sample.per, (params.nperv + 1) * nTarget, 1);
            cnt = cnt + 1;
        end
        
        samplesrecord(tid) = samplesrecord(tid) + 1;
    end
%     dbglist
end

Z.nSamples = cnt - 1;
display(['Average # of trial for one sampling ' num2str(tempcnt/N, '%.02f')]);
samplesrecord
end

function [ar, prior] = getAcceptanceRatio(prevZ, X, sample, tsample, tid, params, prior)
if nargin >= 7
    cachePrior = 1;
else
    cachePrior = 0;
    prior = [];
end
% compute P(Zt|Xt-1)
% copmute P(Xt|Zmt)
lp1 = 0; % log lkhood/prob for prev sample
lp2 = 0; % log lkhood/prob for current sample

if tid == 1
    % motion model.
    for i = 1:size(sample.per, 2)
        % compute observation lkhood
        xid = i;

        ImPred1 = getProjection(sample.per(:, i), sample.cam);
        ImPred2 = getProjection(tsample.per(:, i), tsample.cam);
        
        if xid ~= -1 % is there observation(detection)?
            ImObs = X.obs(:, xid);
            ImObs = [ImObs(1) + ImObs(3) / 2; ImObs(2) + ImObs(4); ImObs(4)];
            temp = ImObs - ImPred1';
            lp1 = lp1 + temp' * (params.Prec1 / ((ImObs(3) / 128) ^ 2)) * temp;
            
            temp = ImObs - ImPred2';
            lp2 = lp2 + temp' * (params.Prec1 / ((ImObs(3) / 128) ^ 2)) * temp;
        end
    end
    % no targets' states are changed/won't contribute on the ratio
    lp1 = exp(-0.5 * lp1) * sum(prod(prior, 1));
    tempcam = repmat(tsample.cam, 1, prevZ.nSamples) - prevZ.cam;
    prior(1, :) = exp(-.5 * diag(tempcam' * params.PrecCam * tempcam))';
    lp2 = exp(-0.5 * lp2) * sum(prod(prior, 1));
else
    tid = tid - 1;
    % compute observation lkhood
    xid = tid;
    ImPred1 = getProjection(sample.per(:, tid), sample.cam);
    ImPred2 = getProjection(tsample.per(:, tid), tsample.cam);
    if xid ~= -1 % is there observation(detection)?
        ImObs = X.obs(:, xid);
        ImObs = [ImObs(1) + ImObs(3) / 2; ImObs(2) + ImObs(4); ImObs(4)];
        temp = ImObs - ImPred1';
        lp1 = lp1 + temp' * (params.Prec1 / ((ImObs(3)/128) ^ 2)) * temp;
        temp = ImObs - ImPred2';
        lp2 = lp2 + temp' * (params.Prec1 / ((ImObs(3)/128) ^ 2)) * temp;
    end
    lp1 = exp(-0.5 * lp1) * sum(prod(prior, 1));
    if tsample.per(6, tid) == 1
        prior(tid + 1, :) = 1/sqrt(2 * pi * 0.1^2) * exp((sample.per(5, tid)-1.7)^2/(-2 * 0.1^2))' * ones(1, prevZ.nSamples);
    else
        prior(tid + 1, :) = 1/3 * ones(1, prevZ.nSamples);
    end
    lp2 = exp(-0.5 * lp2) * sum(prod(prior, 1));
end

ar = lp2 / lp1;

end