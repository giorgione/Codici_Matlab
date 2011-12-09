clear
testlist = [4, 13, 17, 20, 42, 58, 121, 122, 224];

for i = 1:length(testlist)
    [TP, FP, FN, F] = drawTrackingResult(['ResultsFin/test' num2str(testlist(i)) '_-0.20.mat'], 'jpg', '_ant.mat', 1, 2);
    close all
    movie2avi(F(1:2:end), ['ResultVideo/test' num2str(testlist(i)) '_full.avi'], 'FPS', 15);
    clear F;
    
    [TP, FP, FN, F] = drawTrackingResult(['ResultsFin/test' num2str(testlist(i)) '_-0.20_3.mat'], 'jpg', '_ant.mat', 1, 2);
    close all
    movie2avi(F(1:2:end), ['ResultVideo/test' num2str(testlist(i)) '_woG.avi'], 'FPS', 15);
    clear F;
    
    [TP, FP, FN, F] = drawTrackingResult(['ResultsFin/test' num2str(testlist(i)) '_-0.20_5.mat'], 'jpg', '_ant.mat', 1, 2);
    close all
    movie2avi(F(1:2:end), ['ResultVideo/test' num2str(testlist(i)) '_woI.avi'], 'FPS', 15);
    clear F;
end
