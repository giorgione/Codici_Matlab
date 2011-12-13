clear;
addpath(genpath('libs'));

testNums = [400:415];

for i = 1:length(testNums)
    disp(['Process C:/Users/wgchoi/VisionLabWork/HAProject/MATLAB/data/test' num2str(testNums(i)) '/']);
    ProcessKLT(['C:/Users/wgchoi/VisionLabWork/HAProject/MATLAB/data/test' num2str(testNums(i)) '/'], 'jpg', 1300);
%     detectPerson(['data/test' num2str(testNums(i)) '/'], '_det06th.txt', -2);
end



for i = 1:length(testNums)
    showKLT(['C:/Users/wgchoi/VisionLabWork/HAProject/MATLAB/data/test' num2str(testNums(i)) '/'], 'jpg', ['C:/Users/wgchoi/VisionLabWork/HAProject/MATLAB/data/test' num2str(testNums(i)) '/features.txt'])
    pause
end