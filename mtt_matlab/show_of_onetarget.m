function show_of_onetarget(resfile, fstep, kltfile, tid, img_dir, vmode, repeat)

addpath('common\KLTinterface');

[klt_x, klt_y, klt_val] = klt_read_featuretable(kltfile);

load(resfile);
imfiles = dir([img_dir '/*.jpg']);


pid = -1;
for i = 1:length(Tracks)
    if tid == Tracks(i).id
        pid = Tracks(i).tid;
        break;
    end
end

if pid == -1
    disp('ERROR : no matching id!!');
    return;
end

tid = pid;

while(1)
    bblast = [];
    
    for i = 1:fstep:length(imfiles)
        Im = imread([img_dir '/' imfiles(i).name]);

        idx = find(Zs(i).peridx == tid);

        if isempty(idx)
            continue;
        end

        mcam = mean(Zs(i).cam, 2);
        [bb] = getImageProjectionsWRTCam(Zs(i), idx, sparams, mcam);


        figure(1); clf; set(gca,'YDir','reverse')
        title(['Frame #' num2str(i)])
        if vmode == 1
            axis([1 720 1 480])
            pbaspect([3 2 1])

            rectangle('position', bb);

            kidx = find(klt_x(:, i) > bb(1) & klt_x(:, i) < bb(1) + bb(3) & klt_y(:, i) > bb(2) & klt_y(:, i) < bb(2) + bb(4));
            hold on; 
            scatter(klt_x(kidx, i), klt_y(kidx, i), '.'); 
            if ~isempty(bblast)
                for j = kidx'
                    if klt_val(j, i) == 0 && klt_val(j, i - 1) == 0
                        line([klt_x(j, i) klt_x(j, i - fstep)], [klt_y(j, i) klt_y(j, i - fstep)]);
                    end
                end
            end
            hold off;
        else
            axis([1 100 1 200])
            pbaspect([1 2 1])

            kidx = find(klt_x(:, i) > bb(1) & klt_x(:, i) < bb(1) + bb(3) & klt_y(:, i) > bb(2) & klt_y(:, i) < bb(2) + bb(4));

            feats = [klt_x(kidx, i), klt_y(kidx, i)];

            feats = feats - repmat(bb(1:2)', length(kidx), 1);
            feats = feats * 100 / bb(3);

            hold on; 
            scatter(feats(:, 1), feats(:, 2), '.'); 
            
            if ~isempty(bblast)
                diffx = []; diffy = [];
                for j = 1:length(kidx)
                    if klt_val(kidx(j), i) == 0 && klt_val(kidx(j), i - 1) == 0
                        diffx(end + 1) = klt_x(kidx(j), i) - klt_x(kidx(j), i - fstep);
                        diffy(end + 1) = klt_y(kidx(j), i) - klt_y(kidx(j), i - fstep);
                    end
                end
                mmotion = [mean(diffx), mean(diffy)];

%                 mmotion = bb(1:2) - bblast(1:2);
            end
            
            if ~isempty(bblast)
                for j = 1:length(kidx)
                    if klt_val(kidx(j), i) == 0 && klt_val(kidx(j), i - 1) == 0
                        diffx = klt_x(kidx(j), i) - klt_x(kidx(j), i - fstep);
                        diffy = klt_y(kidx(j), i) - klt_y(kidx(j), i - fstep);

                        diffx = diffx * 100 / bb(3); diffy = diffy * 100 / bb(3);
                        line([feats(j, 1) feats(j, 1) - diffx + mmotion(1)], [feats(j, 2) feats(j, 2) - diffy + mmotion(2)]);
                    end
                end
            end
            hold off;
        end
        
        bblast = bb;

        grid on;
        drawnow; % pause(0.3)
    end
    
    if repeat == 0
        break;
    end
end

end

function [bb] = getImageProjectionsWRTCam(Z, id, param, cam)

for sidx = 1:Z.nSamples
    person = Z.per((1 + (id-1) * (param.nperv + 1)):(id * (param.nperv + 1)-1), sidx);
    [dummy, temp] = getProjection(person, cam);
    bb(:, sidx) = temp';
end

end