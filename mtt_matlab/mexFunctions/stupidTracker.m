function stupidTracker(imdir)

imfiles = dir([imdir '*.jpg']);

im1 = imread([imdir imfiles(1).name]);
imshow(im1);
rect = floor(getrect);

patch = im1(rect(2):(rect(2)+rect(4)-1), rect(1):(rect(1)+rect(3)-1), :);
hist = mexComputeAppearanceLkhood(patch, [], 4, .4, 3);

for i = 2:2:length(imfiles)
    im = imread([imdir imfiles(i).name]);
    try
    [histt, rect, sim] = descentSearch(im, hist, rect);
    catch
        break;
    end
    
    if sim < 0.94
        break; 
    end
    
    hist = 0.6 * hist + 0.4 * histt;
    imshow(im);
    rectangle('Position', rect, 'EdgeColor', 'm', 'LineWidth', 2);
    title(['Max sim ' num2str(sim)]);
    drawnow
%     pause;
end

end

function [histout, rectout, sim] = descentSearch(im, hist, rect)
direction = [0, -1, 1; -1, -1, 1; 1, -1, 1; 0, 0, 1; -1, 0, 1; 1, 0, 1; 0, 1, 1; -1, 1, 1; 1, 1, 1];
direction = [direction; 0, -1, .9; -1, -1, .9; 1, -1, .9; 0, 0, .9; -1, 0, .9; 1, 0, .9; 0, 1, .9; -1, 1, .9; 1, 1, .9];
direction = [direction; 0, -1, 1.1; -1, -1, 1.1; 1, -1, 1.1; 0, 0, 1.1; -1, 0, 1.1; 1, 0, 1.1; 0, 1, 1.1; -1, 1, 1.1; 1, 1, 1.1];

temp_hist = [];
sim = -100;

histout = [];
rectout = ones(1, 4);

for i = 1:size(direction, 1)
    temp_rect = floor([rect(1) - direction(i, 1), rect(2) - direction(i, 2), rect(3) * direction(i, 3), rect(4) * direction(i, 3)]);
    patch = im(temp_rect(2):(temp_rect(2)+temp_rect(4)-1), temp_rect(1):(temp_rect(1)+temp_rect(3)-1), :);
    [temp_hist, temp_sim] = mexComputeAppearanceLkhood(patch, hist, 4, .4, 3);
    
%     if temp_sim < 0.85
%         continue;
%     end
    
    if sim < temp_sim
        sim = temp_sim;
        rectout = temp_rect;
        histout = temp_hist;
    end
end

if sum(rectout ~= rect) > 0 && sim ~= -100
    [histout, rectout, sim] = descentSearch(im, hist, rectout);
end

end