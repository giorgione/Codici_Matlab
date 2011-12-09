function [corres] = getCorrespondence2(Z, X, Y, params, sparams, tpan)
if nargin < 6
    tpan = 0;
end

% debug
if Z.nTargets ~= size(Y, 2)
    Z.nTargets
    size(Y, 2)
    error
end

% compute distance betwen existing targets and observation by overlap
distmatrix = zeros(Z.nTargets, length(X.idx));
for j = 1:Z.nTargets
%     [dummy, bbz] = getProjection(zsamp(tidx), csamp);
    [bb] = getImageProjections(Z, j, sparams, 1, tpan);
    bb=mean(bb, 2);
    for k=1:length(X.idx)
        distmatrix(j, k) = distmatrix(j, k) + getOverlap(bb', X.obs(:, k)');
    end
end
% for i = 1:params.notrial
%     sidx = ceil(rand * Z.nSamples);
%     csamp = Z.cam(:, sidx);
%     zsamp = Z.per(:, sidx);
%     
%     for j = 1:Z.nTargets
%         tidx = ((1 + (j - 1) * (params.nperv + 1)):(j * (params.nperv + 1)));
%         
%         [dummy, bbz] = getProjection(zsamp(tidx), csamp);
%         
%         for k=1:length(X.idx)
%             distmatrix(j, k) = distmatrix(j, k) + getOverlap(bbz, X.obs(:, k)') / params.notrial;
%         end
%     end
% end
% compute distance between targets' color hists and detections' color hists
appmatrix = zeros(size(Y, 2), length(X.idx));
for i = 1:size(Y, 2)
    mstrack = Y(:, i);
    if sum(mstrack) ~= 0
        bby = [mstrack(1:2)', mstrack(3)/2, mstrack(3)];
        bby(1:2) = bby(1:2) - bby(3:4) / 2;

        for j=1:length(X.idx)
            appmatrix(i, j) = getOverlap(bby, X.obs(:, j)');
        end
    else
        %
        appmatrix(i, :) = distmatrix(i, :);
    end
end
% for i = 1:size(Y, 2)
%     if 0
%         mstrack = Y(:, i);
%         if sum(mstrack) ~= 0
%             bby = [mstrack(1:2)', mstrack(3)/2, mstrack(3)];
%             bby(1:2) = bby(1:2) - bby(3:4) / 2;
% 
%             for j=1:length(X.idx)
%                 appmatrix(i, j) = getOverlap(bby, X.obs(:, j)');
%             end
%         else
%             %
%             appmatrix(i, :) = distmatrix(i, :);
%         end
%     else
%         for j=1:length(X.idx)
%             try
%                 appmatrix(i, j) = bhattacharyya(Z.model(i).qS(:), X.qS{j}(:));
% %                 appmatrix(i, j) = sum(sum(sum(min(Z.model(i).qS, X.qS{j}))));
%             catch
%                 keyboard
%             end
%         end
%     end
% end
% if sum(sum(appmatrix > 0 & appmatrix < exp(-params.appth))) > 0
%     aaaa= 0;
% end
distmatrix = -log(distmatrix);
appmatrix = -log(appmatrix);

invalididx = (distmatrix > params.ovth & appmatrix > params.appth);
distmatrix(invalididx) = Inf;
appmatrix(invalididx) = Inf;

% distmatrix
% appmatrix

% distmatrix = params.alpha * (b .* distmatrix) + params.beta * (a .* appmatrix);
distmatrix = params.alpha * distmatrix + params.beta * appmatrix;
if sum(sum(isnan(distmatrix))) > 0
    distmatrix(isnan(distmatrix)) = Inf;
end

corrmat = Hungarian(distmatrix);
corres = [];

for i = 1:length(X.idx)
    idx = find(corrmat(:, i) == 1);
    
    if isempty(idx)
        corres(i) = 0;
    else
        corres(i) = idx;
    end
end

end
