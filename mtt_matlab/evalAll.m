testlist = [4, 13, 17, 20, 42, 58, 121, 122, 224];
annstep = [2 2 10 10 10 10 10 10 10];
for i = 1:length(testlist)
    [tTP(i, :), tFP(i, :), tFN(i, :), dTP(i, :), dFP(i, :), dFN(i, :), sl(i)] = evalRecall(['../Dataset/test' num2str(testlist(i)) '/'], 'jpg', ['ResultsFin/test' num2str(testlist(i)) '_'], '.mat', -.5:0.1:0.1, annstep(i));
end

for i = 1:length(testlist)
    [tTP3(i, :), tFP3(i, :), tFN3(i, :)] = evalRecall(['../Dataset/test' num2str(testlist(i)) '/'], 'jpg', ['ResultsFin/test' num2str(testlist(i)) '_'], '_3.mat', -.5:0.1:0.1, annstep(i));
end

for i = 1:length(testlist)
    [tTP5(i, :), tFP5(i, :), tFN5(i, :)] = evalRecall(['../Dataset/test' num2str(testlist(i)) '/'], 'jpg', ['ResultsFin/test' num2str(testlist(i)) '_'], '_5.mat', -.5:0.1:0.1, annstep(i));
end

plot(sum(tFP,1) ./ sum(sl), sum(tTP,1) ./ sum(tTP+tFN, 1), 'ro-', 'LineWidth', 2)
hold on
plot(sum(tFP3,1) ./ sum(sl), sum(tTP3,1) ./ sum(tTP3+tFN3, 1), 'bs-', 'LineWidth', 2)
plot(sum(tFP5,1) ./ sum(sl), sum(tTP5, 1) ./ sum(tTP5 + tFN5, 1), 'cd-', 'LineWidth', 2)
plot(sum(dFP,1) ./ sum(sl), sum(dTP,1) ./ sum(dTP+dFN, 1), 'g^-', 'LineWidth', 2)
hold off
title('Recall/FPPI')
xlabel('FPPI')
ylabel('Recall')
legend({'Full Model', 'No Geometry', 'No Interaction', 'Baseline Detector'})
grid on
axis([0 3 0.6 .9])
