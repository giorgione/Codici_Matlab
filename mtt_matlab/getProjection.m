%function [x] = getProjection(z, cam)
%
% Calculates the projection x of the 3D point z 
% in image Plane according to the camera model.
%
% INPUT:
%                     
%   - cam = [focal, height, xcenter, yhorizon, panning, absolute velocity,];
% 
%   - z = [x, z, vx, vz, heigth] --> 3D Point in Camera System
%
% OUTPUT:
% 
%   - x = [x, y, imhiehgt] --> Image Point
%   
%   - bb = bounding box
%
%  EXTRINSIC CAMERA MODEL:
%
%  Zw: Z in World coordinate system
%  Zc: Z in Camera coordinate system
%
%          Rcw                            Rwc
%           |                              |   -1
%  Zc = R(theta)*Zw + T     <-->  Zw = R(theta)    (Zc-T)    
%                                             |
%                                   (we work with this)
%
%  Rwc: WORLD COORDINATE to CAMERA COORDINATE Rotation Matrix  
%            -1                                    -1
%  Rwc =  Rcw    =   [ cos(cam(5)) ,  -sin(cam(5) ]   
%                    [ sin(cam(5)) ,   cos(cam(5) ] 
%
%  PROJECTION MODEL:
% 
% u = f*x/z + uo
% v = f*y/z + vo   cam(2) = height = y of the classical pinhole model
function [x, bb] = getProjection(z, cam)


R = [cos(cam(5)), sin(cam(5)); 
    -sin(cam(5)), cos(cam(5))];
if(length(cam) == 8)
    z(1:2) = z(1:2) - cam(7:8);
end
%Camera Alignment with World cordinate System
z(1:2) = R * z(1:2);
%z(3) is hz

%Projective Transform to get pixel (u,v)
% u = f*x/z + uo
% v = f*y/z + vo   cam(2) = height = y of the classical pinhole model

%Pixel coordinate unormalized: BOTTOM center point location
x(1) = cam(1) / z(2) * z(1) + cam(3);
x(2) = cam(1) / z(2) * cam(2) + cam(4);
%object height:
x(3) = cam(1) / z(2) * z(5);

% BOUNDING BOX of the Observation
if nargout > 1
    % Upper-Right Point of the bbox 
    %              half-width of the box
    bb(1) = x(1) - x(3) / 4;
    bb(2) = x(2) - x(3);
    
    %width of box
    bb(3) = x(3) / 2;
    %height of box
    bb(4) = x(3);
end

end