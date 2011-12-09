function showDetections(imgdir, imgext, detext, antext)

imfiles = dir([imgdir '/*.' imgext]);
th = 0;
for i = 1:4:length(imfiles)
    imshow([imgdir, imfiles(i).name]);
    fileNameNoExt = imfiles(i).name(1:find(imfiles(i).name=='.',1,'last') - 1);
    if exist([imgdir fileNameNoExt antext])
        load([imgdir fileNameNoExt antext]);
    else
        oneFrameAnnotation = {};
    end

    det = load([imgdir fileNameNoExt detext]);
    det(1, :) = [];
    det(det(:, 6) < th, :) = [];
    
    for j = 1:size(det, 1)
        rectangle('Position', det(j, 1:4), 'LineWidth', 2);
    end
    
    for j = 1:numel(oneFrameAnnotation)
%         rectangle('Position', oneFrameAnnotation{j}.rect, 'LineWidth', 2, 'EdgeColor', 'r');
    end
    
    pause(0.3);
end

end