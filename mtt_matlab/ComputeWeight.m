% version description
% 1. correspondence problem solved by hungarian algorithm
% 2. only one to one matching is allowed
function [ Z, Theta ] = ComputeWeight(Z, Theta, X, Y, corres, params)
% compute weights given correspondence
targetidx = unique(Z.idx);
wtable = zeros(length(Theta.idx), length(Z.idx));

Prec = inv(params.Qdet); % observation noise
Prec2 = inv(params.Qktrack); % ms track noise.
% make a precise equation for the likelihood estimate!!!!!!!!

tic;
temp_idx = [0, targetidx];
corres = corres + 1;


for ts = 1:length(Theta.idx)
    for zs = 1:length(Z.idx)
        % can change it later on.
        xid = find(temp_idx(corres) == Z.idx(zs));
        yid = find(targetidx == Z.idx(zs));
        
        ImPred = getProjection(Z.state(:, zs), Theta.state(:, ts));
        
        tempDistance = 0;
        
        yObs = Y(:, yid);
        if sum(yObs) ~= 0
%             yObs(3:4) = [yObs(3)/2; yObs(3)];
%             yObs(1:2) = yObs(1:2) - yObs(3:4) / 2;
            yObs(2) = yObs(2) + yObs(3) / 2;
            temp = yObs - ImPred';
            tempDistance = temp' * (Prec2 / ((yObs(3) / 128) ^ 2)) * temp;
        end
        
        if ~isempty(xid)
            ImObs = X.obs(:, xid);
            ImObs = [ImObs(1) + ImObs(3) / 2; ImObs(2) + ImObs(4); ImObs(4)];
            temp = ImObs - ImPred';
            
            tempDistance = tempDistance + temp' * (Prec / ((ImObs(3) / 128) ^ 2)) * temp;
        end
        
        if tempDistance == 0
            wtable(ts, zs) = 0;
        else
            wtable(ts, zs) = exp(-1/2 * tempDistance);
        end
    end
end
toc;
% marginalize/normalize over Z and Theta respectively
tempw = sum(wtable, 2);
norm = sum(tempw);

if norm == 0
    Theta.w = (1 / length(Theta.idx) * ones(1, length(Theta.idx)));
else
    Theta.w = (tempw / norm)';
end

for zidx = targetidx(:)'
    id = find(Z.idx == zidx);
    tempw = sum(wtable(:, id), 1);
    norm = sum(tempw);
    
    if norm == 0
        Z.w(id) = (1 / length(id) * ones(1, length(id)));
    else
        Z.w(id) = (tempw / norm)';
    end
end

end
