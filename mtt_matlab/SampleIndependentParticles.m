function [Z] = SampleIndependentParticles(prevZ, params)

targetidx = unique(prevZ.idx);

Z = prevZ;

Z.idx = [];
Z.w = [];
Z.state = [];

avgsamples = 0;
cnt = 1;
for  i = targetidx(:)'
    idx = find(prevZ.idx == i);
    a = cumsum(prevZ.w(idx));
    
    if isfield(params, 'rKernel')
        h = computeEntrophy(prevZ, idx, params.rKernel);
        numPerTarget = ceil(params.asr * h);
        numPerTarget = max(30, numPerTarget);
    else
        numPerTarget = params.numPerTarget;
    end

    avgsamples = avgsamples + numPerTarget;
    tranNoise = mvnrnd(zeros(1, params.lenState), params.R, numPerTarget);
    
    Z.idx = [Z.idx, zeros(1, numPerTarget)];
    Z.w = [Z.w, zeros(1, numPerTarget)];
    Z.state = [Z.state, zeros(params.lenState, numPerTarget)];

    for j = 1:numPerTarget
        samp = sum(a < rand) + 1;
        
        Z.state(:, cnt) = params.A * prevZ.state(:, idx(samp)) + tranNoise(j, :)'; % add noise component.
        Z.w(cnt) = 1;
        Z.idx(cnt) = i;
        
        cnt = cnt + 1;
    end
end
disp(['avg samples : ' num2str(avgsamples/length(targetidx(:)))]);

end


function h = computeEntrophy(z, idx, rKernel)
h = 0;
ndim = size(z.state, 1);
sqKernel = rKernel^2;
for i = idx
    temph = 0;
    for j = idx
        tempd = z.state(:, i) - z.state(:, j);
        temph = temph + z.w(j) * exp(-.5 * dot(tempd, tempd) / sqKernel);
    end
    
    h = h + z.w(i) * log(temph /(sqrt(2 * pi) * rKernel)^ndim);
end
h = -h;
end

% for joint models.
% build a matrix P(interaction_i_j | X_i^{t-1}, X_j^{t-1})
% sample interaction_i_j
% if interaction_i_j == 1 -> consider interaction
% if not ignore it
% do it for every possible pairs.
% after that, sample from P(X_i^t | X_i^{t-1}, X_~i^{t-1})