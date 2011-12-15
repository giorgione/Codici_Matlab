%function [gf] = getGIProj(x, cam)
% Inverse Projection for a point image x=(u,v) on the ground given the Camera
% Status specified by cam.
%
% PARAMETERS INPUT
%
% - x = [u, v] Point on the ground in to IMAGES
% - cam = [focal, height, xcenter, yhorizon];
%
% PARAMETERS OUTPUT
% - gf = [x, z]  3D Point on the ground h=0
function [gf] = getGIProj(x, cam)

gf = zeros(2, size(x, 2));
%INVERSE PROJECTION
gf(2) = cam(1) * cam(2) / (x(2) - cam(4));
gf(1) = (x(1) - cam(3)) * gf(2) / cam(1);

% ROTATION MATRIX
R = [cos(cam(5)), -sin(cam(5)); 
     sin(cam(5)), cos(cam(5))];

%INVERSE ROTATION
gf(1:2) = R * gf(1:2);

%TRANSLATION with Initial Position
if(length(cam) == 8)
    gf(1:2) = gf(1:2) + cam(7:8);
end

end