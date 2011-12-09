function [dTP, dFP, dFN] = computeDetectorRecall(imgdir, detext, annext, annstep, thlist)

detfiles = dir([imgdir '/*' detext]);
dTP = [];
dFP = [];
dFN = [];

for th = thlist
    for i = 1:annstep:length(detfiles)
        fileNameNoExt=detfiles(i).name(1:find(detfiles(i).name=='_',1,'last') - 1);
        if exist([imgdir fileNameNoExt annext])
            load([imgdir fileNameNoExt annext]);
        else
            oneFrameAnnotation = {};
        end
        
        det = load([imgdir detfiles(i).name]);
        det(1, :) = [];
        det(det(:, 6) < th, :) = [];
        
        [TP(i), FP(i), FN(i)] = evalOneFrame(det(:, 1:4)', oneFrameAnnotation);
    end
    
    det(:, 1) = det(:, 1) + det(:, 3) / 8;
    det(:, 3) = det(:, 3) / 4 * 3;
    
    dTP(end + 1) = sum(TP(1:annstep:end));
    dFP(end + 1) = sum(FP(1:annstep:end));
    dFN(end + 1) = sum(FN(1:annstep:end));
end

end


function [TP, FP, FN] = evalOneFrame(bbs, annotations)
removelist = [];
for i = 1:numel(annotations)
    if annotations{i}.rect(4) < 60
        removelist = [removelist, i];
    end
end
annotations(removelist) = [];

if size(bbs, 2) > 0
    % do evaluation!
    olmatrix = [];
    for i = 1:numel(annotations)
        for j = 1:size(bbs, 2);
            olmatrix(i, j) = getOverlap(annotations{i}.rect, bbs(:, j)');
        end
    end
    olmatrix(olmatrix < .55) = 0;
    distmat = -log(olmatrix);

    corrmat = Hungarian(distmat);

    TP = sum(sum(corrmat, 1) == 1);
    FP = sum(sum(corrmat, 1) == 0);
    FN = sum(sum(corrmat, 2) == 0);
else
    TP = 0;
    FP = 0;
    FN = numel(annotations);
end

if (TP + FN) ~= numel(annotations)
    keyboard;
end
end