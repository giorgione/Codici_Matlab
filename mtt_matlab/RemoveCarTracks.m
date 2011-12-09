function [Z, hCTracks, CTracks] = RemoveCarTracks(Z, hCTracks, CTracks, Zidx, sparams)
removelist = [];
for k = 1:length(hCTracks)
    idx = find(hCTracks(k).tid == Z.caridx(Zidx));
    if ~isempty(idx)
        removelist = [removelist, k]; 
    end
end
hCTracks(removelist) = [];

for k = 1:length(CTracks)
    idx = find(CTracks(k).tid == Z.caridx(Zidx));
    if ~isempty(idx)
        CTracks(k).term = i;
    end
end

idx = [];
for j = Zidx
    idx = [idx, ((sparams.ncarv + 1) * (j - 1) + 1):((sparams.ncarv + 1) * j)];
end

Z.car(idx, :) = [];
Z.caridx(Zidx) = [];
Z.tccnt(Zidx) = [];
Z.nCarTargets = Z.nCarTargets  - length(Zidx);


end