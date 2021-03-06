function [z] = getIProjection(x, cam)
% z = [x, z, h, vx, vz]
% x = [u, v, h]
% cam = [focal, height, xcenter, yhorizon];
z = zeros(5, size(x, 2));
z(2) = cam(1) * cam(2) / (x(2) - cam(4));
z(1) = (x(1) - cam(3)) * z(2) / cam(1);
z(3) = 0;
z(4) = 0;
z(5) = x(3)*z(2)/cam(1);

R = [cos(cam(5)), -sin(cam(5)); sin(cam(5)), cos(cam(5))];
z(1:2) = R * z(1:2);

if(length(cam) == 8)
    z(1:2) = z(1:2) + cam(7:8);
end

end