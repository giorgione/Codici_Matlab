function [corres] = getCorrespondence(Z, Theta, X, Im, params)

% weightmask = fspecial('gaussian',[64 64], 21);
weightmask = fspecial('gaussian',[32 64], 10);
weightmask = imresize(weightmask, [128, 64]);
weightmask = weightmask(:) * 1000; % prevent error due to small value max(w) * 1000 = 0.4...
edges = 32 * (0:8);

% construct histograms for observations
for i=1:length(X.idx)
    patch = getImgPatch(Im, X.obs(:, i)');
%     patch = Im(X.obs(2, i):(X.obs(2, i) + X.obs(4, i) - 1), X.obs(1, i):(X.obs(1, i) + X.obs(3, i) - 1), :);
    patch = imresize(patch, [128 64]);
    patch = patch(:); patch = [patch(1:128*64), patch((128*64+1):(128*64*2)), patch((128*64*2+1):(128*64*3))];
    a = patch < 0 | patch > 255; patch(a) = 0;
    h = histc_nD(patch, edges, weightmask); h=h(:);
    X.h(:, i) = h;
end

% compute distance betwen existing targets and observation by overlap
% dist = -log(bboverlap)
% targetidx = unique(Z.idx);
targetidx = unique(Z.idx);

distmatrix = zeros(length(targetidx), length(X.idx));
a = cumsum(Theta.w);
for i = 1:length(targetidx)
    zidx = find(Z.idx == targetidx(i));
    b = cumsum(Z.w(zidx));
    
    for k = 1:params.notrial
        csamp = sum(a < rand * a(end)) + 1;
        zsamp = sum(b < rand * b(end)) + 1;

        [dummy, bbz] = getProjection(Z.state(:, zidx(zsamp)), Theta.state(:, csamp));
        
        for j=1:length(X.idx)
            distmatrix(i, j) = distmatrix(i, j) + getOverlap(bbz, X.obs(:, j)') / params.notrial;
        end
    end
end
%%%%
% for i = 1:length(Z.pest)
%     [dummy, bbz] = getProjection(Z.pest(i).loc, Theta.pest.loc);
%     for j=1:length(X.idx)
%         distmatrix(i, j) = getOverlap(bbz, X.obs(:, j)');
%     end
% end
%%%%

distmatrix = -log(distmatrix);
distmatrix(distmatrix > params.ovth) = Inf;

% compute distance between targets' color hists and detections' color hists
appmatrix = zeros(length(targetidx), length(X.idx));
for i = 1:length(targetidx)
    for j=1:length(X.idx)
        appmatrix(i, j) = bhattacharyya(Z.model(:, i), X.h(:, j));
    end
end
% apply gating on appearance histogram
appmatrix(appmatrix > params.appth) = Inf;

distmatrix
appmatrix

% before adding up we have to check hundered percent certain cases....

a = distmatrix > 0.01;
b = appmatrix > 0.01;

% distmatrix = params.alpha * (b .* distmatrix) + params.beta * (a .* appmatrix);
distmatrix = params.alpha * distmatrix + params.beta * appmatrix;
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

% update app model.
for i = 1:length(targetidx)
    a = find(corres == i);
    if isempty(a)
        continue;
    end
    Z.model(:, i) = (params.adf) * Z.model(:, i) + (1 - params.adf) * X.h(:, a);
end
t = 1;

end

%function ol = getOverlap(A, B)
% Calculates the OVERLAP between the bounding box A and B
% PARAMETERS INPUT:
%
% A: Detection Vector [x, y, width, height]
% B: Detection Vector [x, y, width, height]
%
% ol:  OVERLAP between the bounding box of Detections A and B
%
%       o-------------------o A1
%       |o-------------o B1 |
%       ||             |    |  
%       ||             |    |     
%       ||             |    |
%       o|-------------|----o
%       A|             |    
%        o-------------o    
%        B                          
%
function ol = getOverlap(A, B)
  
%calculate upper right corner fo the bbox of Detections
A1 = A + [0, 0, A(1:2) - 1];
B1 = B + [0, 0, B(1:2) - 1];

%Calculate overlap box
obox = [max(A1(1), B1(1)), max(A1(2), B1(2)), min(A1(3), B1(3)), min(A1(4), B1(4))];
TotArea=(A(3) * A(4) + B(3) * B(4));
ol = 2 * max(0, obox(3)-obox(1)+1) * max(0, obox(4)-obox(2)+1) /TotArea ;
end