function [x] = getProjectionVector(z, cam)

% cam = [focal, height, xcenter, yhorizon]';
% x = [x, y, imhiehgt]'
% z = [x, z, vx, vz, heihgt]'

x(1, :) = cam(1, :) ./ z(2, :) .* z(1, :) + cam(3, :);
x(2, :) = cam(1, :) ./ z(2, :) .* cam(2, :) + cam(4, :);
x(3, :) = cam(1, :) ./ z(2, :) .* z(5, :);

end