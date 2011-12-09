function [TP, FP, FN] = showOneFrame(Im, frame, Z, Tracks, CTracks, oneFrameAnnotation, sparams, KLT, tKLT, det, Xc, Y, mode)
targetidx = unique(Z.peridx);
caridx = unique(Z.caridx);
figure(1);
imshow(Im);
figure(2);
clf;
xbound = 15;
ccc = 'rgbck';
lst = {'-', '-', ':', '-.'};

isz = size(Im);

if frame > 1
    figure(1)
    validklts = [];
    drawidx = find(Z.gfcnt > 1);
    for i = 1:length(drawidx)
        validklts = [validklts, find(tKLT.idx == Z.gfidx(drawidx(i)))];
    end
    hold on
    scatter(tKLT.x, tKLT.y, 'filled', 'g');
    scatter(tKLT.x(validklts), tKLT.y(validklts), 'filled', 'r');
    hold off
    
    %%%%%%%%%% draw ground plane.
    if 0
        mcam = mean(Z.cam, 2);

        hold on
        [xx, zz] = getVisibleRegion(mcam);
        for i = 1:size(xx, 1)
            [x1] = getGProj([xx(i, 1), xx(i, 2)]', mcam);
            [x2] = getGProj([xx(i, 1), xx(i, 3)]', mcam);

            plot(linspace(x1(1), x2(1), 10), linspace(x1(2), x2(2), 10), 'm:');
        end
        for i = 1:size(zz, 1)
            [x1] = getGProj([zz(i, 2), zz(i, 1)]', mcam);
            [x2] = getGProj([zz(i, 3), zz(i, 1)]', mcam);

            plot(linspace(x1(1), x2(1), 10), linspace(x1(2), x2(2), 10), 'm:');
        end
        hold off 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if 1
    figure(2)
    hold on
    scatter(mean(Z.gfeat((drawidx-1) * 3 + 1, :), 2), mean(Z.gfeat((drawidx-1) * 3 + 2, :), 2), 'rx')
    tic;
    hold off
    end
end
if 1
if sparams.useInteraction
    figure(2);
    hold on
    tempbeta = mean(Z.beta, 4);
    for j = 1:size(Z.beta, 1)
        trackid = 0;
        for l = 1:length(Tracks)
            if (Tracks(l).term == -1) && (Z.peridx(j) == Tracks(l).tid)
                trackid = l;
                break;
            end
        end
        if trackid == 0, continue; end

        for k = (j+1):size(Z.beta, 1)
            trackid = 0;
            for l = 1:length(Tracks)
                if (Tracks(l).term == -1) && (Z.peridx(k) == Tracks(l).tid)
                    trackid = l;
                    break;
                end
            end
            if trackid == 0, continue; end

            [dummy, id] = max(tempbeta(j, k, :));
            if id == 1
                continue;
            end

            person1 = mean(Z.per((1 + (j-1) * (sparams.nperv + 1)):(j * (sparams.nperv + 1)), :), 2);
            person2 = mean(Z.per((1 + (k-1) * (sparams.nperv + 1)):(k * (sparams.nperv + 1)), :), 2);

            if id == 2
                plot([person1(1) person2(1)], [person1(2) person2(2)], 'm', 'LineWidth', 2);
            end
        end
    end
    hold off
end
end

if mode == 0
    evalmatch = zeros(10, length(oneFrameAnnotation));
    for j = 1:sparams.nretry
        sidx = j;
        figure(1);
        rectangle('Position', [1, Z.cam(4, sidx), isz(2), 1], 'LineStyle', ':', 'EdgeColor', 'r');
        if 1
        figure(2);
        hold on;
        va = [xbound; 2 * Z.cam(1, sidx) / isz(2) * xbound];
        R = [cos(Z.cam(5, sidx)), -sin(Z.cam(5, sidx)); sin(Z.cam(5, sidx)), cos(Z.cam(5, sidx))];
        va = R * va + [Z.cam(7:8, sidx)];
        plot([Z.cam(7, sidx) va(1)], [Z.cam(8, sidx) va(2)], 'r:', 'LineWidth', 2);

        va = [-xbound; 2 * Z.cam(1, sidx) / isz(2) * xbound];
        va = R * va + [Z.cam(7:8, sidx)];
        plot([Z.cam(7, sidx) va(1)], [Z.cam(8, sidx) va(2)], 'r:', 'LineWidth', 2);

        %%% draw camera!
        hold off;
        end
        tempcnt = 0;
        for k = targetidx(:)'
            trackid = 0;
            for l = 1:length(Tracks)
                if (Tracks(l).term == -1) && (k == Tracks(l).tid)
                    trackid = l;
                    break;
                end
            end
            if trackid == 0, continue; end

            if sum(Tracks(trackid).occlusion == frame) > 0
%                 continue;
            end
            
            zidx = find(Z.peridx == Tracks(trackid).tid);
            if Z.per((sparams.nperv + 1) * zidx, sidx) <= 0.1
                continue;
            end

            person = Z.per((1 + (zidx-1) * (sparams.nperv + 1)):(zidx * (sparams.nperv + 1)), sidx);
            [dummy, bb] = getProjection(person, Z.cam(:, sidx));
            
            figure(1);
            rectangle('Position', bb, 'LineWidth', 2, 'EdgeColor', ccc(mod(trackid, 5)+1), 'LineStyle',  lst{mod(ceil(trackid / 5), 4) + 1});
            hold on;
            if j == 1
                text(bb(1) + 3, bb(2) + bb(4) - 15, num2str(trackid), 'Fontsize', 25, 'Color', ccc(mod(trackid, 5)+1), 'FontWeight', 'bold');
            end
            quiver(bb(1) + bb(3)/2, bb(2) + bb(4)/2, person(3) * 30, -person(4) * 30, [ccc(mod(trackid, 5)+1)], 'LineWidth', 2);
            hold off;

            if 0
            figure(2);
            hold on
            scatter(person(1), person(2), 'filled', [ccc(mod(trackid, 5)+1)]);
            if j == 1
                text(person(1), person(2) + 1, num2str(trackid), 'Fontsize', 13, 'Color', ccc(mod(trackid, 5)+1), 'FontWeight', 'bold');
            end
            quiver(person(1), person(2), person(3), person(4), ccc(mod(trackid, 5) + 1), 'LineWidth', 2);
            hold off
            end
            % evaluation routine!
            for l = 1:length(oneFrameAnnotation)
                evalmatch(tempcnt + 1, l) = evalmatch(tempcnt+ 1, l) + getOverlap(oneFrameAnnotation{l}.rect, bb);
            end
            tempcnt = tempcnt + 1;
        end
        
        for k = caridx(:)'
            trackid = 0;
            for l = 1:length(CTracks)
                if (CTracks(l).term == -1) && (k == CTracks(l).tid)
                    trackid = l;
                    break;
                end
            end
            if trackid == 0, continue; end

            zidx = find(Z.caridx == CTracks(trackid).tid);
            if Z.car((sparams.ncarv + 1) * zidx, sidx) <= 0.1
                continue;
            end

            car = Z.car((1 + (zidx-1) * (sparams.ncarv + 1)):(zidx * (sparams.ncarv + 1)), sidx);
            [dummy, bb] = getProjectionCar(car, Z.cam(:, sidx));
            
            figure(1);
            rectangle('Position', bb, 'LineWidth', 1, 'EdgeColor', ccc(mod(trackid, 5)+1), 'LineStyle',  lst{mod(ceil(trackid / 5), 4) + 1});
            hold on;
            if j == 1
                text(bb(1) + 3, bb(2) + bb(4) - 15, num2str(trackid), 'Fontsize', 25, 'Color', ccc(mod(trackid, 5)+1), 'FontWeight', 'bold');
            end
            hold off;
        end
    end
else
    
end

% evaluation routine!
if tempcnt > 0
    TP = sum(sum(evalmatch(1:tempcnt, :) / sparams.nretry > .5, 2) > 0);
    FP = sum(sum(evalmatch(1:tempcnt, :) / sparams.nretry > .5, 2) == 0);
    FN = sum(sum(evalmatch(1:tempcnt, :) / sparams.nretry > .5, 1) == 0);
else
    TP = 0;
    FP = 0;
    FN = length(oneFrameAnnotation);
end

figure(1);
for j = 1:size(det, 1)
    rectangle('Position', det(j, 1:4), 'LineStyle', ':', 'LineWidth', 1, 'EdgeColor', 'k');
end

for j = 1:size(Xc.obs, 2)
    rectangle('Position', Xc.obs(:, j), 'LineStyle', ':', 'LineWidth', 2, 'EdgeColor', 'k');
end

for j = 1:size(Y, 2)
    if sum(Y(:, j)) == 0
        continue;
    end
    bbox = [Y(1, j) - Y(3, j) / 4, Y(2, j) - Y(3, j) / 2, Y(3, j) / 2, Y(3, j)];
    rectangle('Position', bbox, 'LineStyle', ':', 'LineWidth', 1, 'EdgeColor', 'w');
end

title(['frame' num2str(frame)]);

drawnow;
if 1
figure(2)
% ref = [floor((mcam(7) + 2.5)/5) * 5; floor((mcam(8) + 2.5)/5) * 5]
% axis([ref(1) - xbound, ref(1) + xbound,  ref(2) - 1, ref(2) + 2*xbound*3/4]);
% axis([-xbound xbound 0 2*xbound*3/4]);
grid on;
pbaspect([1 2 1])
drawnow;
end
end


function [xx, zz] = getVisibleRegion(mcam)

xref = getIProjection([320; 480; 1], mcam);
xcenter = floor(xref(1) / 5 + .5) * 5;
xx(:, 1) = ((xcenter - 15):5:(xcenter + 15))';
xx(:, 2) = xref(2);
if cos(mcam(5)) > 0
    xx(:, 3) = 1000;
else
    xx(:, 3) = -1000;
end
zcenter = floor(xref(2) / 5 + .5) * 5;


end