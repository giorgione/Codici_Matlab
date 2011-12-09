function [gf] = getGIProj(x, cam)
% gf = [x, z]
% x = [u, v]
% cam = [focal, height, xcenter, yhorizon];
gf = zeros(2, size(x, 2));
gf(2) = cam(1) * cam(2) / (x(2) - cam(4));
gf(1) = (x(1) - cam(3)) * gf(2) / cam(1);

R = [cos(cam(5)), -sin(cam(5)); sin(cam(5)), cos(cam(5))];
gf(1:2) = R * gf(1:2);


if(length(cam) == 8)
    gf(1:2) = gf(1:2) + cam(7:8);
end

end