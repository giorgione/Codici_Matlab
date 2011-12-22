%function [Z, hTracks, Tracks] = RemoveTracks(Z, hTracks, Tracks, Zidx, sparams)
% Rimuovi da hTrack le track specificate da Zidx

function [Z, hTracks, Tracks] = RemoveTracks(Z, hTracks, Tracks, Zidx, sparams)
removelist = [];
for k = 1:length(hTracks)
    idx = find(hTracks(k).tid == Z.peridx(Zidx));
    if ~isempty(idx)
        removelist = [removelist, k]; 
    end
end
hTracks(removelist) = [];

for k = 1:length(Tracks)
    idx = find(Tracks(k).tid == Z.peridx(Zidx));
    if ~isempty(idx)
        Tracks(k).term = idx; %era i
    end
end

idx = [];
for j = Zidx
    idx = [idx, ((sparams.nperv + 1) * (j - 1) + 1):((sparams.nperv + 1) * j)];
end
Z.per(idx, :) = [];
Z.model(Zidx) = [];
Z.peridx(Zidx) = [];
Z.tcnt(Zidx) = [];
Z.nTargets = Z.nTargets - length(Zidx);

%%%%% remove interaction too!!!!!
Z.beta(:, Zidx, :, :) = [];
Z.beta(Zidx, :, :, :) = [];

end