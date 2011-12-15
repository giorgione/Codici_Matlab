%function Z = FilterOutGFeats(Z, list, sparams)
% Removes GROUND FEATURES specified in list from Z
function Z = FilterOutGFeats(Z, list, sparams)
    
if Z.nFeats > 0
        %Generate the index list to be removed
        removeidx = [];
        for j = list
            removeidx = [removeidx, [(1:sparams.ngfeat) + (j-1)*sparams.ngfeat]];
        end
        %remove 
        Z.gfeat(removeidx, :) = [];
       
        %remove the ground features VALUES in list from Z.gfV
        Z.gfV(list, :) = [];
        %remove the ground features IDX in list Z.gfidx
        Z.gfidx(list) = [];
        %remove the ground features VALUES in list Z.gfcnt
        Z.gfcnt(list) = [];
        
        Z.nFeats = Z.nFeats - length(list);
    end
end
