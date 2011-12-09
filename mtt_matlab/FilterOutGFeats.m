function Z = FilterOutGFeats(Z, list, sparams)
if Z.nFeats > 0
    removeidx = [];
    for j = list
        removeidx = [removeidx, [(1:sparams.ngfeat) + (j-1)*sparams.ngfeat]];
    end
    Z.gfeat(removeidx, :) = [];
    Z.gfV(list, :) = [];
    Z.gfidx(list) = [];
    Z.gfcnt(list) = [];
    Z.nFeats = Z.nFeats - length(list);
end
end
