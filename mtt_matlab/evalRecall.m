function [tTP, tFP, tFN, dTP, dFP, dFN, seqLen] = evalRecall(filedir, imgext, prefix, postfix, list, fstep)
tTP = [];
tFP = [];
tFN = [];
dTP = [];
dFP = [];
dFN = [];
for i = list
    disp(['Processing ' [prefix num2str(i, '%.02f') postfix]]);
    [TP, FP, FN] = drawTrackingResult([prefix num2str(i, '%.02f') postfix], imgext, '_ant.mat', 0, fstep);
    tTP(end+1) = sum(TP(1:fstep:end));
    tFP(end+1) = sum(FP(1:fstep:end));
    tFN(end+1) = sum(FN(1:fstep:end));
end

if nargin >= 4
    seqLen = length(FN(1:fstep:end));
    [dTP, dFP, dFN] = computeDetectorRecall(filedir, '_det06th.txt', '_ant.mat', fstep, list);
end

end