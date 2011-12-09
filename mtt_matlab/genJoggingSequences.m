clear

% track targets
testNums = 400%:415;
% testNums = 400;
ihor = [];

dirprefix = 'C:\Users\wgchoi\VisionLabWork\HAProject\MATLAB\data\test';
imgext = 'jpg';

for i = testNums
    dirname = [dirprefix num2str(i) '\'];
    files = dir([dirname '*.' imgext]);
    imshow([dirname files(1).name]);
    [x, y] = getpts;
    ihor(end + 1) = y(1);
end

a = 1;

for i = 1:length(testNums)
    dirname = [dirprefix num2str(testNums(i)) '\'];
    dirname
    ihor(i)
    TrackTop(dirname, 11, ihor(i), ['ss_test' num2str(testNums(i))], -0.1);
end