% function [corres] = getCorrespondence2(Z, X, Y, params, sparams, tpan)
%
% Generate the first part of the AFFINITY MATRIX
%
%
function [corres] = getCorrespondence2(Z, X, Y, params, sparams, tpan)
%Set the panning angle if unspecified
if nargin < 6
    tpan = 0;
end

% debug
if Z.nTargets ~= size(Y, 2)
    Z.nTargets
    size(Y, 2)
    error
end

% DISTANCE MATRIX: distance betwen existing targets and observation 
% by overlap
distmatrix = zeros(Z.nTargets, length(X.idx));
% For each OBSERVATION j of (DETECTED TARGET) estimated in the previous frame:
%
%                                         Distr. over the target j
%                               ^          |
%1) Project Z    in the image:  X     = E[ X      | Z       , THETA    ] 
%            (j,t)               (j,t)      (j,t)    (j,t-1)      (t-1) 
%                                  |                    |             |
%                    estimate of TARGET         Estimated 3D status  camera 
%                    j at time t in the         value of TARGET j    conf
%                    image                      a time t-1           t-1
%
%2) Obtain an estimate of the bounding box ,                                       
%                                                 ^
%3) Compute the INTERSECTION BETWEN bounding box  X     U    X
%                                                  (j,t)     (j,t)
%4) Computea affinity matrix
for j = 1:Z.nTargets
    %BBOX of the ESTIMATED TARGET in the previous frame
    [bb] = getImageProjections(Z, j, sparams, 1, tpan);
    
    % Expected Values for the PREDICTION of TARGET j
    bb=mean(bb, 2);
    
    %OVERLAP RATIO A(j,k)
    for k=1:length(X.idx)
        distmatrix(j, k) = distmatrix(j, k) + getOverlap(bb', X.obs(:, k)');
    end
end
 

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
 

%Compute LOGARITM
distmatrix = -log(distmatrix);
appmatrix = -log(appmatrix);

%Invalid Association must be set to Inf to garantue that Hungarian
%alghoritm will not use them for the final solution
invalididx = (distmatrix > params.ovth & appmatrix > params.appth);
distmatrix(invalididx) = Inf;
appmatrix(invalididx) = Inf;
 
%Calculate the final AFFINITY MATRIX: weighted sum distmatrix and appmatrix
distmatrix = params.alpha * distmatrix + params.beta * appmatrix;
if sum(sum(isnan(distmatrix))) > 0
    distmatrix(isnan(distmatrix)) = Inf;
end

%Apply Hungarian Alghoritm to the Association Problem:
% binary matrix
corrmat = Hungarian(distmatrix);
corres = [];

%Generate the Vector of Correspondeces :
%               0   --  if i Object hasn't correspondence
% corres(i) = 
%              idx  of the Correspondent Object to i
for i = 1:length(X.idx)
    %idx of the Correspondent Dete
    idx = find(corrmat(:, i) == 1);
    
    if isempty(idx)
        corres(i) = 0; %No correspondence for i Observation
    else
        corres(i) = idx; %Save the idx
    end
end

end
