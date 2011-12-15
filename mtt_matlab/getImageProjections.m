% function [bb] = getImageProjections(Z, id, param, pred, campan)
%  Return the bounding box of the Detected Persons for each sampled Camera STATUS
% 
function [bb] = getImageProjections(Z, id, param, pred, campan)
if nargin < 4
    %First frame processing
    pred = 0;
    campan = 0;
else if nargin < 5
    %Other frames and camera pannning angle not specified
     campan = 0;
end
%Consider each GENERATED SAMPLE for the PERSON VARIABLE
% sidx: sample index
for sidx = 1:Z.nSamples
    ind1=(1 + (id-1) * (param.nperv + 1));
    ind2=id * (param.nperv + 1)-1; 
    
    %Extract the person  data from Z
    person = Z.per((ind1:ind2), sidx);
    
    %Generate the STATE VALUES for CAMERA VARIABLE cam
    
    if pred == 1
        % When the previous frame is the First frame of the video sequence 
        % Z is [] so we need to :
        % 1) extract the CAMERA STATUS from the Evolution of CAMERA
        % 2) extract the Persono from the Evolution
        % Evolution of Person
        person = param.Aper * person;
        %Camera Evolution
        cam = getNextCamera(Z.cam(:, sidx), param);
        %Update CAMERA PAN
        cam(5) = cam(5) + campan;
    else
        %Other frames: 
        %Extract the Last one camera STATUS from Z
        cam = Z.cam(:, sidx);
    end
    
    %Project the 3D-PERSON on the 2D image according to the camera Status
    [dummy, temp] = getProjection(person, cam);
    %Person BBOX 
    bb(:, sidx) = temp';
end
end