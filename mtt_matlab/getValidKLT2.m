function idx = getValidKLT2(KLT, i, Z, imwidth, sparams, pred)
if nargin < 6
    pred = 1;
end

%%%%%%% we may want to reduce number of KLT features...
mcam = mean(Z.cam, 2);
idx = find(KLT.x(:, i) > sparams.KLTmargin & KLT.y(:, i) > mcam(4) + 30 & KLT.x(:, i) < imwidth - sparams.KLTmargin);
for j = 1:length(Z.peridx)
    [bb] = getImageProjections(Z, j, sparams, pred);
    idx = FeatsNotInBB(bb, KLT, idx, i);
end

end

function idx = FeatsNotInBB(bb, KLT, idx, i)
margin = 20;
inidx = find(KLT.x(idx, i) > (bb(1) - margin) & KLT.x(idx, i) < (bb(1) + bb(3) + margin) & KLT.y(idx, i) > (bb(2) - margin) & KLT.y(idx, i) < (bb(2) + bb(4) + margin));
idx(inidx) = [];
end