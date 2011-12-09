function idx = getValidKLT(KLT, i, Z, imwidth, nperv)
i = i - 1;
if i == 0
    idx = [];
    return;
end
mcam = mean(Z.cam, 2);
idx = find(KLT.val(:, i + 1) == 0 & KLT.y(:, i) > mcam(4) - 10 & KLT.x(:, i) > 60 & KLT.x(:, i) < imwidth - 60);
for trial = 1:5
    sidx = ceil(rand * Z.nSamples);
    
    for zidx = 1:length(Z.peridx)
        if Z.per((nperv + 1) * zidx, sidx) == 0
            continue;
        end

        person = Z.per((1 + (zidx-1) * (nperv + 1)):(zidx * (nperv + 1)), sidx);
        [dummy, bb] = getProjection(person, Z.cam(:, sidx));
        
        rectangle('position', bb);
        
        idx = FeatsNotInBB(bb, KLT, idx, i);
    end
end

end

function idx = FeatsNotInBB(bb, KLT, idx, i)
margin = 10;
inidx = find(KLT.x(idx, i) > (bb(1) - margin) & KLT.x(idx, i) < (bb(1) + bb(3) + margin) & KLT.y(idx, i) > (bb(2) - margin) & KLT.y(idx, i) < (bb(2) + bb(4) + margin));
idx(inidx) = [];
end