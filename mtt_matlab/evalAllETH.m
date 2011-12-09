testlist = [2 3];
annstep = [4 1];

for i = 1:length(testlist)
    [tTP(i, :), tFP(i, :), tFN(i, :), dTP(i, :), dFP(i, :), dFN(i, :), sl(i)] = evalRecall(['../Dataset/seq0' num2str(testlist(i)) '/'], 'png', ['ResultsFin/seq0' num2str(testlist(i)) '_'], '.mat', -.7:0.1:0.2, annstep(i));
end


for i = 1:length(testlist)
    [tTP2(i, :), tFP2(i, :), tFN2(i, :)] = evalRecall(['../Dataset/seq0' num2str(testlist(i)) '/'], 'png', ['ResultsFin/seq0' num2str(testlist(i)) '_'], '_4.mat', [-.5 -.2 .1], annstep(i));
end


for i = 1:length(testlist)
    [tTP3(i, :), tFP3(i, :), tFN3(i, :)] = evalRecall(['../Dataset/seq0' num2str(testlist(i)) '/'], 'png', ['ResultsFin/seq0' num2str(testlist(i)) '_'], '_6.mat', [-.5 -.2 .1], annstep(i));
end