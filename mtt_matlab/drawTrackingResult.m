function [TP, FP, FN, F] = drawTrackingResult(resfile, imgext, annext, bshow, fstep, pref, img_dir)
if  nargin < 4
    bshow = 1;    img_dir = [];    fstep = 1; pref = [];
elseif nargin < 5
    img_dir = []; fstep = 1; pref = [];
elseif nargin < 6
    img_dir = []; pref = [];
elseif nargin < 7
    img_dir = [];
end

if nargout > 3
    record = 1;
else
    record = 0;    F = [];
end
load(resfile);
% 'Zs', 'Tracks', 'TP', 'FP', 'FN', 'sparams',  'imgdir',  'framerate', 'KLTused'
if isempty(img_dir)
    img_dir = imgdir;
end
imfiles = dir([img_dir '/*.' imgext]);
cmap = colormap;

for i = 1:length(Tracks)
    Tracks(i).occlusion = [];
    for frame = 1:length(imfiles)
        if (sum(Tracks(i).det > frame - 7 & Tracks(i).det <= frame) >= 1)
            Tracks(i).initFrame = frame;
            break;
        end
    end
end

TP = zeros(1, length(imfiles)); 
FP = zeros(1, length(imfiles));
FN = zeros(1, length(imfiles));

xbound = 15;

camtraj = [];
for i = 1:fstep:min(length(imfiles), length(Zs))
    camtraj = [camtraj, Zs(i).cam(7:8)];
end

%%%%%%%%%%%%%
OvelayTrak = 0;

for i = 1:fstep:min(length(imfiles), length(Zs))
    Im = imread([img_dir '/' imfiles(i).name]);
    [Tracks] = inferOcclusion(Tracks, Zs(i), i, sparams);
    fileNameNoExt=imfiles(i).name(1:find(imfiles(i).name=='.',1,'last') - 1);
    if exist([img_dir fileNameNoExt annext])
        load([img_dir fileNameNoExt annext]);
%         truecnt(i) = length(oneFrameAnnotation);
    else
%         truecnt(i) = 0; % truecnt(i - fstep);
        oneFrameAnnotation = {};
    end
    
    if bshow == 1
        figure(1);
        imshow(Im);
    end
    isz = size(Im);
    
    [T1, T2, T3, wtraj] = drawTargetsOnImagePlane(Zs, Tracks, i, cmap, sparams, oneFrameAnnotation, bshow, fstep, isz);
    
    [wtraj_c] = drawCarTargetsOnImagePlane(Zs, CTracks, i, cmap, sparams, bshow, fstep, isz);
    TP(i)=T1; FP(i)=T2; FN(i)=T3;
    if bshow == 1
        if exist('KLTused')
            drawKLTOnImagePlane(Zs, KLTused, i, sparams);
        end
        if 0 
            if exist([img_dir fileNameNoExt '_det06th.txt'])
                det = load([img_dir fileNameNoExt '_det06th.txt']);
                det(1, :) = [];
                det(det(:, 6) < 0, :) = [];
            end

            for j = 1:size(det, 1)
                rectangle('position', det(j, 1:4), 'linewidth', 2);
            end
        end
        text(10, 25, ['Frame ' num2str(i)], 'Fontsize', 20, 'Color', 'k', 'FontWeight', 'bold');

        figure(2);
        va = [xbound; 2 * Zs(i).cam(1) / isz(2) * xbound];
        R = [cos(Zs(i).cam(5)), -sin(Zs(i).cam(5)); sin(Zs(i).cam(5)), cos(Zs(i).cam(5))];
        va = R * va + [Zs(i).cam(7:8)];
        plot([Zs(i).cam(7) va(1)], [Zs(i).cam(8) va(2)], 'r:', 'LineWidth', 5);
        hold on;
        va = [-xbound; 2 * Zs(i).cam(1) / isz(2) * xbound];
        va = R * va + [Zs(i).cam(7:8)];
        plot([Zs(i).cam(7) va(1)], [Zs(i).cam(8) va(2)], 'r:', 'LineWidth', 5);
        plot(camtraj(1, 1:fstep:floor(i/fstep)), camtraj(2, 1:fstep:floor(i/fstep)), 'r:', 'LineWidth', 5);


        for j = 1:length(wtraj)
            cid = getColIdx(wtraj(j).id);
            plot(wtraj(j).traj(1, :), wtraj(j).traj(2, :), 'Color', cmap(cid, :), 'LineWidth', 3);
            plot(wtraj(j).pos(1), wtraj(j).pos(2), '.', 'MarkerEdgeColor', cmap(cid, :), 'MarkerFaceColor', cmap(cid, :), 'MarkerSize', 25);
        end
        
        for j = 1:length(wtraj_c)
            cid = getColIdx(wtraj_c(j).id);
            plot(wtraj_c(j).traj(1, :), wtraj_c(j).traj(2, :), 'Color', cmap(cid, :), 'LineWidth', 3);
            plot(wtraj_c(j).pos(1), wtraj_c(j).pos(2), 'd', 'MarkerEdgeColor', cmap(cid, :), 'MarkerFaceColor', cmap(cid, :), 'MarkerSize', 15);
        end

        gfeats = reshape(Zs(i).gfeat, 3, length(Zs(i).gfeat)/3);

        scatter(gfeats(1, :), gfeats(2, :), 100, 'r.');
%         scatter(Zs(i).);
        %%% draw camera!


        if sparams.useInteraction
            tempbeta = mean(Zs(i).beta, 4);
            for j = 1:size(Zs(i).beta, 1)
                trackid = 0;
                for l = 1:length(wtraj)
                    if Zs(i).peridx(j) == wtraj(l).tid
                        trackid = l;
                        break;
                    end
                end
                if trackid == 0, continue; end
                for k = (j+1):size(Zs(i).beta, 1)
                    trackid = 0;
                    for l = 1:length(wtraj)
                        if Zs(i).peridx(k) == wtraj(l).tid
                            trackid = l;
                            break;
                        end
                    end
                    if trackid == 0, continue; end

                    [dummy, id] = max(tempbeta(j, k, :));
                    if id == 1
                        continue;
                    end

                    person1 = mean(Zs(i).per((1 + (j-1) * (sparams.nperv + 1)):(j * (sparams.nperv + 1)), :), 2);
                    person2 = mean(Zs(i).per((1 + (k-1) * (sparams.nperv + 1)):(k * (sparams.nperv + 1)), :), 2);

                    if id == 2
                        plot([person1(1) person2(1)], [person1(2) person2(2)], 'm', 'LineWidth', 2);
                    end
                end
            end
        end

        grid on;
        hold off;
        title(['Frame ' num2str(i)]);
%             axis([-15 15 0 25])
        axis([Zs(i).cam(7)-20 Zs(i).cam(7)+20 Zs(i).cam(8)-20 Zs(i).cam(8)+20])

        drawnow;
        
        if record == 1
            figure(1);
            F(i) = getframe;
            F(i).cdata = imresize(F(i).cdata, [240 320]);
            figure(2);
            F2 = getframe;
            F2.cdata = imresize(F2.cdata, [240 240]);
            
            tcdata(:,:,1) = [F(i).cdata(:,:,1), F2.cdata(:,:,1)];
            tcdata(:,:,2) = [F(i).cdata(:,:,2), F2.cdata(:,:,2)];
            tcdata(:,:,3) = [F(i).cdata(:,:,3), F2.cdata(:,:,3)];
            
            F(i).cdata = tcdata;
        end
    end
end

end

function [wtraj] = drawCarTargetsOnImagePlane(Zs, CTracks, frame, cmap, sparams, bshow, fstep, isz)

wtraj = [];
mcam = mean(Zs(frame).cam, 2);

bbs = [];
bdraw = zeros(1, length(CTracks));

for i = 1:length(CTracks)
    if (CTracks(i).det(1) > frame)
        continue;
    end
    
    if (CTracks(i).det(length(CTracks(i).det)) <= frame - 6)
        continue;
    end

    zidx = find(Zs(frame).caridx == CTracks(i).tid);
    if length(zidx) == 0
        continue;
    end

    bdraw(i) = 1;
    
    tempwt.pos = Zs(frame).car((7 * zidx - 6):(7 * zidx - 5));
    
    tempwt.id = CTracks(i).id;
    tempwt.tid = CTracks(i).tid;
    tempwt.traj = zeros(2, 0);
    
    [bb] = getImageProjectionsCars(Zs(frame), zidx, sparams);
    bb = mean(bb, 2);
    
    if (bb(1) + bb(3) * 1 / 3 < 0) || ...
        (bb(1) + bb(3) * 2 / 3 >= isz(2)) || ...
        (bb(2) + bb(4) * 2 / 3 >= isz(1))
        continue;
    end
            
    bbs = [bbs, bb];
    
    if bshow == 1
        cid = getColIdx(CTracks(i).id);
        rectangle('Position', bb, 'Linewidth', 5, 'EdgeColor', cmap(cid, :));
        text(bb(1) + 3, bb(2) + bb(4) - 15, num2str(CTracks(i).id), 'Fontsize', 9, 'Color', cmap(cid, :));

        traj = [bb(1) + bb(3) / 2; bb(2) + bb(4)];
        for j = 2:fstep:100
            if sum(CTracks(i).det <= frame - j) == 0
                break;
            end

            zidx = find(Zs(frame-j).caridx == CTracks(i).tid);

            tempwt.traj = [tempwt.traj, Zs(frame-j).car((7 * zidx - 6):(7 * zidx - 5))];
%             bb = getImageProjectionsWRTCam(Zs(frame-j), zidx, sparams, mcam);
%             bb = mean(bb, 2);
% 
%             if (bb(2) + bb(4)) < 150
%                 break;
%             end
%             traj(:, end+1) = [bb(1) + bb(3) / 2; bb(2) + bb(4)];
        end

        if size(traj, 2) > 0
%             plot(traj(1, :), traj(2, :), 'LineWidth', 2, 'Color', cmap(cid, :));
        end
    end
    
    wtraj = [wtraj, tempwt];
end

end

function [TP, FP, FN, wtraj] = drawTargetsOnImagePlane(Zs, Tracks, frame, cmap, sparams, annotations, bshow, fstep, isz)
if nargin < 6
    annotations = {};
    eval = 0;
else
    eval = 1;
    % filter out annotations smaller than 60 pixels following the
    % convention of ETH
    removelist = [];
    for i = 1:numel(annotations)
        if annotations{i}.rect(4) < 60
            removelist = [removelist, i];
        end
    end
    annotations(removelist) = [];
end
wtraj = [];

mcam = mean(Zs(frame).cam, 2);
if bshow == 1
    hold on
    plot([1 720], [mcam(4) mcam(4)], 'r:', 'LineWidth', 2);
end

bbs = [];
bdraw = zeros(1, length(Tracks));
for i = 1:length(Tracks)
    % validcheck
%     Tracks(i)
% try
    if (Tracks(i).initFrame > frame)
        continue;
    end
    
    if (sum(Tracks(i).occlusion == frame) > 0)
        continue;
    end
    
    if (Tracks(i).det(length(Tracks(i).det)) <= frame - 6)
        continue;
    end
% catch ee
%     keyboard
% end
%     disp(['TID : ' num2str(i) ' at ' num2str(frame)]);
    zidx = find(Zs(frame).peridx == Tracks(i).tid);
    if length(zidx) == 0
        continue;
    end

    bdraw(i) = 1;
    
    tempwt.pos = Zs(frame).per((6 * zidx -5):(6 * zidx -4));
    tempwt.id = Tracks(i).id;
    tempwt.tid = Tracks(i).tid;
    tempwt.traj = zeros(2, 0);
    
    [bb] = getImageProjections(Zs(frame), zidx, sparams);
    bb(1) = bb(1) + bb(3) / 8;
    bb(3) = bb(3) / 4 * 3;
    
    bb = mean(bb, 2);
    
    if (bb(1) + bb(3) * 1 / 3 < 0) || ...
        (bb(1) + bb(3) * 2 / 3 >= isz(2)) || ...
        (bb(2) + bb(4) * 2 / 3 >= isz(1))
        continue;
    end
            
    bbs = [bbs, bb];
    
    if bshow == 1
        cid = getColIdx(Tracks(i).id);
%         try
            rectangle('Position', bb, 'Linewidth', 5, 'EdgeColor', cmap(cid, :));
            text(bb(1) + 3, bb(2) + bb(4) - 15, num2str(Tracks(i).id), 'Fontsize', 9, 'Color', cmap(cid, :));


            traj = [bb(1) + bb(3) / 2; bb(2) + bb(4)];
            for j = 2:2:100
                if sum(Tracks(i).det <= frame - j) == 0
                    break;
                end

                zidx = find(Zs(frame-j).peridx == Tracks(i).tid);
                
                tempwt.traj = [tempwt.traj, Zs(frame-j).per((6 * zidx -5):(6 * zidx -4))];
                bb = getImageProjectionsWRTCam(Zs(frame-j), zidx, sparams, mcam);
                bb = mean(bb, 2);

                if (bb(2) + bb(4)) < 150
                    break;
                end
                traj(:, end+1) = [bb(1) + bb(3) / 2; bb(2) + bb(4)];
            end

            if size(traj, 2) > 0
                plot(traj(1, :), traj(2, :), 'LineWidth', 2, 'Color', cmap(cid, :));
            end
%         catch
%         end
    end
    
    wtraj = [wtraj, tempwt];
end
% 
% match = zeros(1, Zs(frame).nTargets);
% 
% for i = 1:Zs(frame).nTargets
%     for k =1:length(Tracks)
%         Zs(frame).peridx()
%     end
% end
% 
% if bshow == 1
%     for i = 1:Zs(frame).nTargets
%         if 
%         for j = (i+1):Zs(frame).nTargets
%             
%         end
%     end
% end
if bshow == 1
    hold off
end


TP = 0;
FP = 0;
FN = 0;

if eval == 1
    if size(bbs, 2) > 0
        % do evaluation!
        olmatrix = [];
        for i = 1:numel(annotations)
            for j = 1:size(bbs, 2);
                olmatrix(i, j) = getOverlap(annotations{i}.rect, bbs(:, j)');
            end
        end
        olmatrix(olmatrix < .55) = 0;
        distmat = -log(olmatrix);

        corrmat = Hungarian(distmat);

        TP = sum(sum(corrmat, 1) == 1);
        FP = sum(sum(corrmat, 1) == 0);
        FN = sum(sum(corrmat, 2) == 0);
    else
        TP = 0;
        FP = 0;
        FN = numel(annotations);
    end
    
    if (TP + FN) ~= numel(annotations)
        keyboard;
    end
end
end

function drawKLTOnImagePlane(Zs, KLTused, frame, sparams)
hold on
scatter(KLTused(frame).tKLT.x, KLTused(frame).tKLT.y, 'filled', 'r');
if length(KLTused(frame).invalid) > 0
    idx = [];
    for i = KLTused(frame).invalid
        idx = [idx, find(KLTused(frame).tKLT.idx == i)];
    end
    scatter(KLTused(frame).tKLT.x(idx), KLTused(frame).tKLT.y(idx), 'filled', 'g');
end

hold off
end

function cidx = getColIdx(tid)
cidx = mod(tid * 7, 64) + 1;
end


function [bb] = getImageProjectionsWRTCam(Z, id, param, cam)

for sidx = 1:Z.nSamples
    person = Z.per((1 + (id-1) * (param.nperv + 1)):(id * (param.nperv + 1)-1), sidx);
    [dummy, temp] = getProjection(person, cam);
    bb(:, sidx) = temp';
end

end