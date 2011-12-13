function getSiftMatchForVideo(imdir, imext, outfile)

addpath(genpath('sift-0.9.19.tar'));

imfiles = dir([imdir '*.' imext]);

mkdir([imdir 'siftmatches/']);

for i = 997:2:(length(imfiles)-2)
    if(exist([imdir 'siftmatches/matches_' num2str(i) '_' num2str(i+2) '.mat']))
        continue;
    end
    
    disp(['processing frame ' num2str(i) ' and ' num2str(i + 2)]);
    while(~exist([imdir 'tempsift/sift_' num2str(i) '.mat']))
        pause(1);
    end
    
    berror = 1;
    while(berror)
        try
            load([imdir 'tempsift/sift_' num2str(i) '.mat'], 'frames', 'descr');%, 'gss', 'dogss');
            berror = 0;
        catch
            berror = 1;
        end
    end
    
    frames1 = frames;
    descr1 = descr;
    
    while(~exist([imdir 'tempsift/sift_' num2str(i + 2) '.mat']))
        pause(1);
    end
    
    berror = 1;
    while(berror)
        try
            load([imdir 'tempsift/sift_' num2str(i+2) '.mat'], 'frames', 'descr');%, 'gss', 'dogss');
            berror = 0;
        catch
            berror = 1;
        end
    end
    
    frames2 = frames;
    descr2 = descr;
    
    matches=siftmatch(descr1, descr2);
    
    save([imdir 'siftmatches/matches_' num2str(i) '_' num2str(i+2) '.mat'], 'matches');
    
    imshow([imdir imfiles(i+2).name]);
    title(['match between ' num2str(i) ' and ' num2str(i + 2)])
    
    hold on;
    for j = 1:size(matches, 2)
        line([frames1(1, matches(1,j)), frames2(1, matches(2,j))], [frames1(2, matches(1,j)), frames2(2, matches(2,j))]);
    end
    hold off;
    
    disp(['found ' num2str(size(matches, 2)) ' number of matches!']);
    
    drawnow
    pause(1);
end

end