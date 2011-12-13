function visualizeKT2(KT2, CT, img_dir, showImage, outputvideo)

global params;

if nargin < 5
    recordVid = 0;
else
    recordVid = 1;
    prevSize = [];
    
end

fl = dir([img_dir '*.jpg']);

Im = imread([img_dir fl(1).name]);
params = [400, 1.5, size(Im, 2) / 2, size(Im, 1) / 2];

% posestr={'N' 'R' 'RF' 'F' 'LF' 'L' 'LB' 'B' 'RB'};
posestr={'N' 'R' 'F' 'L' 'B'};
ccode = 'kbmgcr';
nccode = length(ccode);

WIDTH = 64; HEIGHT = 128;
BBSIZE = [WIDTH;HEIGHT];

minx = -25;
maxx = 25;
maxz = 50;

figHandle = figure;
% figHandle2 = figure;
for i = 1:length(fl)
    if showImage == 1
        figure(figHandle);

        fprintf('showing %dth frame\n', i);
        Im = imread([img_dir fl(i).name]);
        imshow(Im);

        params(2) = CT.X(1, i);
        params(4) = CT.X(2, i);
        
        hold on
        plot([1, size(Im, 2)], [params(4), params(4)], 'r:', 'linewidth', 3);
        hold off
        for j=1:numel(KT2)  %detectionList{i}
            if i < KT2(j).detection_list(1) || i > KT2(j).detection_list(end) %not yet detected, or after termination
                continue;
            end

            color = ccode(mod(j, nccode) + 1);
            
            id = find(KT2(j).detection_list == i);
            
            if ~isempty(id)
                bb = KT2(j).Dets(1:4, id);
                rectangle('Position', bb', 'LineWidth', 1, 'EdgeColor', color, 'LineStyle', ':');
            end
            
            S = KT2(j).State(:, i - KT2(j).detection_list(1) + 1);
%             S = KT2(j).X3(:, i - KT2(j).detection_list(1) + 1);
            if KT2(j).type == 1
                X3 = S(1:3, 1);
                Xp = hfunction(X3, params);
                scale = Xp(3) / HEIGHT;
                
                if (scale <= 0)
                    continue;
                end

                if sum(isnan(Xp)) > 0
                    continue;
                end

                bb = [Xp(1) - WIDTH * scale / 2; Xp(2) - HEIGHT * scale ; BBSIZE * scale];

                if S(3, 1) > 2.3 || S(3, 1) < 1.2
                    LS = ':';
                else
                    LS = '-';
                end
                rectangle('Position', bb', 'LineWidth', 3, 'EdgeColor', color, 'LineStyle', LS);
            else
                X3 = S(1:4, 1);
                Xp = hfunction2(X3, params);
                
                if (sum(Xp < 0) > 0)
                    continue;
                end

                if sum(isnan(Xp)) > 0
                    continue;
                end

                bb = [Xp(1)-Xp(3)/2;Xp(2)-Xp(4); Xp(3:4)];

                if S(4, 1) > 5 || S(4, 1) < 0.5
                    LS = ':';
                else
                    LS = '-';
                end
                rectangle('Position', bb', 'LineWidth', 3, 'EdgeColor', color, 'LineStyle', LS);
            end
        end
        drawnow;
        
        F(i) = getframe;
                
        figure(figHandle);
        clf;
        for j=1:numel(KT2)  %detectionList{i}
            if i< KT2(j).detection_list(1) || i>KT2(j).detection_list(end) %not yet detected, or after termination
                continue;
            end

            color = ccode(mod(j, nccode) + 1);
            
            S = KT2(j).State(1:2, (i - KT2(j).detection_list(1) + 1):-1:1);
            hold on
            if KT2(j).type == 1
                plot(S(1, 1),S(2, 1), 'o', 'Color', color);
            else
                plot(S(1, 1),S(2, 1), '<', 'Color', color);
            end
            plot(S(1, :),S(2, :), ':', 'Color', color);
            text(S(1) + .5, S(2) - 0.1, num2str(j), 'Color', color, 'FontSize', 12, 'FontWeight', 'bold');
            hold off
        end
        axis([minx, maxx, 0, maxz]);
        drawnow;
        
        tF2 = getframe;
        
        stf1 = size(F(i).cdata);
%         F(i).cdata = F(i).cdata(1:480, 1:720, :);
        stf2 = size(tF2.cdata);
        F(i).cdata(1:stf1(1), stf1(2)+1:stf1(2)*2, 1) = [tF2.cdata(:,:,1), zeros(stf2(1), stf1(2) - stf2(2)); zeros(stf1(1) - stf2(1), stf1(2))];
        F(i).cdata(1:stf1(1), stf1(2)+1:stf1(2)*2, 2) = [tF2.cdata(:,:,2), zeros(stf2(1), stf1(2) - stf2(2)); zeros(stf1(1) - stf2(1), stf1(2))];
        F(i).cdata(1:stf1(1), stf1(2)+1:stf1(2)*2, 3) = [tF2.cdata(:,:,3), zeros(stf2(1), stf1(2) - stf2(2)); zeros(stf1(1) - stf2(1), stf1(2))];
    elseif showImage == 2
    elseif showImage == 3
        figure(figHandle);

        fprintf('showing %dth frame\n', i);
        Im = imread([img_dir fl(i).name]);
        imshow(Im);

        params(2) = CT.X(1, i);
        params(4) = CT.X(2, i);
        
        hold on
        plot([1, size(Im, 2)], [params(4), params(4)], 'r:', 'linewidth', 3);
        hold off
        for j=1:numel(KT2)  %detectionList{i}
            if i < KT2(j).detection_list(1) || i > KT2(j).detection_list(end) %not yet detected, or after termination
                continue;
            end

            color = ccode(mod(j, nccode) + 1);
            
            id = find(KT2(j).detection_list == i);
            
            if ~isempty(id)
                bb = KT2(j).Dets(1:4, id);
                rectangle('Position', bb', 'LineWidth', 1, 'EdgeColor', color, 'LineStyle', ':');
            end
            
            S = KT2(j).State(:, i - KT2(j).detection_list(1) + 1);
%             S = KT2(j).X3(:, i - KT2(j).detection_list(1) + 1);
            if KT2(j).type == 1
                X3 = S(1:3, 1);
                Xp = hfunction(X3, params);
                scale = Xp(3) / HEIGHT;
                
                if (scale <= 0)
                    continue;
                end

                if sum(isnan(Xp)) > 0
                    continue;
                end

                bb = [Xp(1) - WIDTH * scale / 2; Xp(2) - HEIGHT * scale ; BBSIZE * scale];

                if S(3, 1) > 2.3 || S(3, 1) < 1.2
                    LS = ':';
                else
                    LS = '-';
                end
                rectangle('Position', bb', 'LineWidth', 3, 'EdgeColor', color, 'LineStyle', LS);
            else
                X3 = S(1:4, 1);
                Xp = hfunction2(X3, params);
                
                if (sum(Xp < 0) > 0)
                    continue;
                end

                if sum(isnan(Xp)) > 0
                    continue;
                end

                bb = [Xp(1)-Xp(3)/2;Xp(2)-Xp(4); Xp(3:4)];

                if S(4, 1) > 5 || S(4, 1) < 0.5
                    LS = ':';
                else
                    LS = '-';
                end
                rectangle('Position', bb', 'LineWidth', 3, 'EdgeColor', color, 'LineStyle', LS);
            end
        end
        drawnow;
    elseif showImage == 4
    else
        figure(figHandle);
        clf;
        for j=1:numel(KT2)  %detectionList{i}
            if i< KT2(j).detection_list(1) || i>KT2(j).detection_list(end) %not yet detected, or after termination
                continue;
            end

            color = ccode(mod(j, nccode) + 1);
            
            S = KT2(j).State(1:2, (i - KT2(j).detection_list(1) + 1):-1:1);
            hold on
            if KT2(j).type == 1
                plot(S(1, 1),S(2, 1), 'o', 'Color', color);
            else
                plot(S(1, 1),S(2, 1), '<', 'Color', color);
            end
            plot(S(1, :),S(2, :), ':', 'Color', color);
            text(S(1) + .5, S(2) - 0.1, num2str(j), 'Color', color, 'FontSize', 12, 'FontWeight', 'bold');
            hold off
        end
        axis([minx, maxx, 0, maxz]);
        drawnow;
    end
    
%     if recordVid == 1
%         F(i) = getframe;
%         if length(prevSize) > 0
%             F(i).cdata = F(i).cdata(1:prevSize(1), 1:prevSize(2), :);
%         else
%             prevSize = size(F(i).cdata);
%         end
%     end
end

if recordVid == 1
     movie2avi(F, outputvideo);
end
end