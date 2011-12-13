function getSiftForVideo(imdir, imext, outfile)

addpath(genpath('sift-0.9.19.tar'));

imfiles = dir([imdir '*.' imext]);

mkdir([imdir 'tempsift']);

for i = 1001:2:length(imfiles)
    try
        Im = imread([imdir imfiles(i).name]);
    catch
        imfiles(i+1:end) = [];
        break;
    end
    
    Im = imresize(Im, .5);
    
    Im = rgb2gray(Im);
    Im = Im-min(Im(:));
    Im = Im/max(Im(:));
    
    [frames,descr,gss,dogss] = sift(Im, 'Verbosity', 1);
    
    frames(1:3, :) = 2 .* frames(1:3, :);
    
    save([imdir 'tempsift/sift_' num2str(i) '.mat'], 'frames', 'descr');%, 'gss', 'dogss');
end


end

