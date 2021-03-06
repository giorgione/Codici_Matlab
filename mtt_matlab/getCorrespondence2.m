% function [corres] = getCorrespondence2(Z, X, Y, params, sparams, tpan)
%
% Solve the  corrispondences between DETECTIONs in current time (in X) 
% and PREDICTED 3D target status in the previous frame (in Z.per) projected
% in the image
% How to:
% For each OBSERVATION j in X (DETECTED TARGET) and estimated TARGET Z 
% in the previous frame:
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
%4) Compute distance matrix Proj(Z) - X
%
% For each TRACKED POSITION j in Y (MATCHED TARGET) and estimated TARGET Z 
% in the previous frame:
%
%1) Compute AFFINITY matrix Y - X
%
% Appy Hungarian to AFFINITY+DISTANCE and get Z-X MATCHED
%
% INPUT:
% 
%  - Z: TRACKED TARGETs in the previous FRAME
%
%  - X: Detected TARGETs (CANDIDATES) in the currente FRAME( Old-New cars/person)
%
%  - Y: Output from the Meanshift Tracker 
% OUTPUT:
%
% -corres: correspondece Vector containing the idx of the Detected Object X
%          to which correspond the 
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

% DISTANCE MATRIX distance by bbox-overlap betwen:
% Z: existing targets 
% X: Detected observation by DETECTORS
distmatrix = zeros(Z.nTargets, length(X.idx));
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

% APPEARENCE MATRIX between 
% - Y: MEAN SHIFT OUTPUT 
%
% - X: DETECTIONS
%
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

% Righe:   Z-Y Index
% Colonne: Detections Index
%Calculate the final ASSOCIATION MATRIX: weighted sum distmatrix and appmatrix
distmatrix = params.alpha * distmatrix + params.beta * appmatrix;
if sum(sum(isnan(distmatrix))) > 0
    distmatrix(isnan(distmatrix)) = Inf;
end

%Apply Hungarian Alghoritm to the Association Problem:
% binary matrix
corrmat = Hungarian(distmatrix);
corres = [];

%Generate the Vector of Correspondeces :
%               0   --  if i DETECTION i hasn't correspondence in Z
% corres(i) = 
%              idx  of the TARGET in Z Correspondent to the i-DETECTIONs in X  
for i = 1:length(X.idx)
    %idx of the TARGET from the previous time in CORRISPONDENCE with
    %DETECTIONS in X
    idx = find(corrmat(:, i) == 1);
    
    if isempty(idx)
        corres(i) = 0; %No correspondence for i Observation
    else
        corres(i) = idx; %Save the idx
    end
end

end
