% function idx = getValidKLT2(KLT, i, Z, imwidth, sparams, pred)
% Search VALID KLT features in the i-Frame and MATCH with KLT of the previous
% frames in Z.
% Validate KLT features in the CURRENT FRAME by FILTERING OUT  KLT points 
% KLT.x(:, i) and KLT.y(:, i) that are:
% 
% 1) OUTSIDE the sub-region R of the image: 
%    * sparams.KLTmargin < KLT.x(:, i) < imwidth - sparams.KLTmargin
%    * KLT.y(:, i) is > mcam(4) <--> horizon line + 30 pixel
%
% 2) INSIDE the BOUNDING-BOX of DETECTED PEOPLE
%
% PARAMETERS INPUT:
%
% - KLT: KLT FEATURES for all the frames--> KLT in the current FRAME
% - i: current FRAME Index
% - Z: PREVIOUS Z configuration --> KLT of the PREVIOUS FRAME 
% - imwidth: image widht
% - sparams: Model Data
% - pred: Previous frame
%
% PARAMETERS OUTPUT:
%
% - idx: INDEX of VALID KLT Features

function idx = getValidKLT2(KLT, i, Z, imwidth, sparams, pred)
    %When the first frame is processed
    if nargin < 6
        pred = 1;
    end

    %%%%%%% we may want to reduce number of KLT features..
    %
    %MONTE CARLO APPROXIMATION: Generate the CAMERA STAUS from Samples
    mcam = mean(Z.cam, 2);
    %Search for VALID KLT Features in i-frame. THEY ARE VALID if:
    %
    % a)  sparams.KLTmargin < x < imwidth - sparams.KLTmargin
    % b)  y is > mcam(4) --> horizon line + 30 pixel
    idx = find(KLT.x(:, i) > sparams.KLTmargin & KLT.y(:, i) > mcam(4) + 30 & KLT.x(:, i) < imwidth - sparams.KLTmargin);

    % If are Detected Person try to filter out KLT features inside their 
    % BOUNDING BOX 
    for j = 1:length(Z.peridx)
        %Extract bbox of Detected Person: KLT in that box ARE NOT VALID
        [bb] = getImageProjections(Z, j, sparams, pred);
        %Get Index of KLT OUTSIDE of that box
        idx = FeatsNotInBB(bb, KLT, idx, i);
    end

end

function idx = FeatsNotInBB(bb, KLT, idx, i)
    margin = 20;
    %get index of KLT-features inside the BBOX  bb
    inidx = find(KLT.x(idx, i) > (bb(1) - margin) & KLT.x(idx, i) < (bb(1) + bb(3) + margin) & KLT.y(idx, i) > (bb(2) - margin) & KLT.y(idx, i) < (bb(2) + bb(4) + margin));
    %delete inidx from idx
    idx(inidx) = [];
end