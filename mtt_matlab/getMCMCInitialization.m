function sample = getMCMCInitialization(X, mcam, params)

sample.cam = mcam;
sample.cam(4) = params.mpcam;
sample.per = [];
for i = 1:size(X.obs, 2)
    % compute inverse projection
    % randomly sample the humanness by looking at the heights...
    x = X.obs(:, i);
    x = [x(1) + x(3) / 2, x(2) + x(4), x(4)];
    sample.per = [sample.per, [getIProjection(x', sample.cam); 1]]; % rand < X.pobj(i)]];
end

M = 1 + size(X.obs, 2);
N = 1000;
horpert = 2;

for i = 1:N
    tsample = sample;
    id = ceil(rand * M);
    if id == 1
        tsample.cam(4) = tsample.cam(4) + randn * horpert;
        if (params.useCamAnnealing == 1)
            horpert = horpert * params.AnnealingConst;
        end
    else
        if rand > 0.9
            tsample.per(end, id-1) = ~tsample.per(end, id-1);
        end
        tsample.per(1:5, id - 1) = tsample.per(1:5, id - 1) + mvnrnd(zeros(1, 5), diag([0.1 0.1 1e-10 1e-10 0.01].^2))';
    end
    
    ar = getIAcceptanceRatio(X, sample, tsample, id, params);
    if rand < ar
        sample = tsample;
    end
end

end

function [ar] = getIAcceptanceRatio(X, sample, tsample, id, params)
lp1 = 0; 
lp2 = 0;

if id == 1
    for i = 1:size(sample.per, 2)
        % compute observation lkhood
        lp1 = lp1 + getObservationLKhood(sample.cam, sample.per(:, i), X.obs(:, i), params);
        lp2 = lp2 + getObservationLKhood(tsample.cam, tsample.per(:, i), X.obs(:, i), params);
    end
    
    lp1 = exp(lp1-lp2);
    lp2 = 1;
else
    i = id - 1;
    lp1 = lp1 + getObservationLKhood(sample.cam, sample.per(:, i), X.obs(:, i), params);
	lp2 = lp2 + getObservationLKhood(tsample.cam, tsample.per(:, i), X.obs(:, i), params);
    
    lp1 = exp(lp1-lp2);
    lp2 = 1;
    if sample.per(6, i) == 1
       lp1 = lp1 * normpdf(sample.per(5, i), 1.7, 0.1);
    else
       lp1 = lp1 * params.noHprob;
    end
    
    if tsample.per(6, i) == 1
       lp2 = lp2 * normpdf(tsample.per(5, i), 1.7, 0.1);
    else
       lp2 = lp2 * params.noHprob;
    end
end

ar = lp2 / (lp1 + 1e-200);
% if abs(lp2/lp1 - ar) > 0.000001
%     keyboard
% end

end

function lkhood = getObservationLKhood(cam, state, Xobs, params)
lkhood = 0;
ImPred = getProjection(state, cam);

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

end