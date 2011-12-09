function [feats] = sift_read_featuretable(imgdir, max_frames)

siftdir = [imgdir 'tempsift/'];
matchdir = [imgdir 'siftmatches/'];

siftfiles = dir([siftdir '/*.mat']);

tempfeats = {};
lastidx = [];

siftfiles(max_frames:end) = [];

for i = 1:2:(length(siftfiles)-2)
    load([siftdir 'sift_' num2str(i) '.mat'], 'frames');
    frames1 = frames;
    load([siftdir 'sift_' num2str(i+2) '.mat'], 'frames');
    frames2 = frames;
    
    load([matchdir 'matches_' num2str(i) '_' num2str(i+2) '.mat'], 'matches');

    % prefiltering
    temp=[];
    for j = 1:size(matches, 2)
        d = sum((frames1(1:2, matches(1, j)) - frames2(1:2, matches(2, j))).^2)^.5;
        if d > 150
            temp(end+1) = j;
        end
    end
    matches(:, temp) = [];
    
    temp=[];
    for j = 1:size(matches, 2)
        idx = find(matches(1, :) == matches(1, j));
        if length(idx) == 1
            continue;
        end
        
        minidx = -1;
        mindist = 10000;
        for k = idx(:)'
            d = sum((frames1(1:2, matches(1, j)) - frames2(1:2, matches(2, k))).^2)^.5;
            if d < mindist
                minidx = k;
                mindist = d;
            end
        end
        
        temp = [temp, setdiff(idx(:)', minidx)];
    end
    
    matches(:, temp) = [];
    temp=[];
    for j = 1:size(matches, 2)
        idx = find(matches(2, :) == matches(2, j));
        if length(idx) == 1
            continue;
        end
        
        minidx = -1;
        mindist = 10000;
        for k = idx(:)'
            d = sum((frames2(1:2, matches(2, j)) - frames1(1:2, matches(1, k))).^2)^.5;
            if d < mindist
                minidx = k;
                mindist = d;
            end
        end
        
        temp = [temp, setdiff(idx(:)', minidx)];
    end
    matches(:, temp) = [];
    
    curidx = -1*ones(length(lastidx), 1);
    for j = 1:size(matches, 2)
%         line([frames1(1, matches(1,j)), frames2(1, matches(2,j))], [frames1(2, matches(1,j)), frames2(2, matches(2,j))]);
        idx = find(lastidx == matches(1, j));
        
        if isempty(idx)
            % 
            curidx(end + 1) = matches(2, j);
            feats.loc = sparse(2,1);
            feats.val = frames1(4, matches(1, j));
            
            feats.loc(:, i) = frames1(1:2, matches(1, j));
            feats.loc(:, i+2) = frames2(1:2, matches(2, j));
            
            tempfeats(end+1) = {feats};
        else
            curidx(idx) = matches(2, j);
            try
            tempfeats{idx}.loc(:, i) = frames1(1:2, matches(1, j));
            
            tempfeats{idx}.loc(:, i+2) = frames2(1:2, matches(2, j));
            
            catch
                keyboard
            end
        end
    end
    
    lastidx = curidx;
end

sz = 0;
for i = 1:length(tempfeats)
    if sz < size(tempfeats{i}.loc, 2)
        sz = size(tempfeats{i}.loc, 2);
    end
end

feats.x = sparse(length(tempfeats),sz);
feats.y = sparse(length(tempfeats),sz);
feats.val = sparse(length(tempfeats),1);
feats.vidx = 1:length(tempfeats);

for i = 1:length(tempfeats)
    feats.x(i, 1:length(tempfeats{i}.loc(1, :))) = tempfeats{i}.loc(1, :);
    feats.y(i, 1:length(tempfeats{i}.loc(1, :))) = tempfeats{i}.loc(2, :);
    feats.val(i) = tempfeats{i}.val;
end
% 
% imfiles = dir([imgdir '*.bmp']);
% for i = 3:2:length(imfiles)
%     subplot(121);
%     imshow([imgdir imfiles(i).name]);
%     
%     hold on;
%     for j = 1:length(tempfeats)
%         if size(tempfeats{j}.loc, 2) >= i && (sum(tempfeats{j}.loc(:, i - 2)) ~= 0 && sum(tempfeats{j}.loc(:, i)) ~= 0)
%             line([tempfeats{j}.loc(1, i-2) tempfeats{j}.loc(1, i)], [tempfeats{j}.loc(2, i-2) tempfeats{j}.loc(2, i)]);
%         end
%     end
%     hold off;
%     
%     subplot(122);
%     imshow([imgdir imfiles(i).name]);
%     
%     load([siftdir 'sift_' num2str(i-2) '.mat'], 'frames');
%     frames1 = frames;
%     load([siftdir 'sift_' num2str(i) '.mat'], 'frames');
%     frames2 = frames;
%     
%     load([matchdir 'matches_' num2str(i-2) '_' num2str(i) '.mat'], 'matches');
%     
%     hold on;
%     for j = 1:size(matches, 2)
%         line([frames1(1, matches(1,j)), frames2(1, matches(2,j))], [frames1(2, matches(1,j)), frames2(2, matches(2,j))]);
%     end
%     hold off;
%     
%     pause;
% end

end