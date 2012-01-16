%function [corres] = getCorrespondenceCars(Z, Xc, params, sparams, tpan)
%
% Solve the association problem between the CAR observation Xc and the CAR
% Tracked in the previous frame (Non integrato con il Meanshift per aggiungere 
% la Matrice di APPARENZA)

function [corres] = getCorrespondenceCars(Z, Xc, params, sparams, tpan)
if nargin < 6
    tpan = 0;
end

% compute distance betwen existing targets and observation by overlap
distmatrix = zeros(Z.nCarTargets, length(Xc.idx));
for j = 1:Z.nCarTargets
    [bb] = getImageProjectionsCars(Z, j, sparams, 1, tpan);
    
    bb=mean(bb, 2);
    
    for k=1:length(Xc.idx)
        distmatrix(j, k) = getOverlap(bb', Xc.obs(:, k)');
    end
end

distmatrix = -log(distmatrix);
invalididx = (distmatrix > params.covth);
distmatrix(invalididx) = Inf;

if sum(sum(isnan(distmatrix))) > 0
    distmatrix(isnan(distmatrix)) = Inf;
end

corrmat = Hungarian(distmatrix);
corres = [];

for i = 1:length(Xc.idx)
    idx = find(corrmat(:, i) == 1);
    
    if isempty(idx)
        corres(i) = 0;
    else
        corres(i) = idx;
    end
end

end