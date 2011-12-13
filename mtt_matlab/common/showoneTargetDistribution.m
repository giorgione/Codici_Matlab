function map = showoneTargetDistribution(m, W, V, idx)

if length(idx) ~= 2
    error;
end

xmin = min(m(idx(1), :)) - sqrt(V{1}(idx(1), idx(1))) * 4;
xmax = max(m(idx(1), :)) + sqrt(V{1}(idx(1), idx(1))) * 4;
ymin = min(m(idx(2), :)) - sqrt(V{1}(idx(2), idx(2))) * 4;
ymax = max(m(idx(2), :)) + sqrt(V{1}(idx(2), idx(2))) * 4;

x = xmin:(xmax - xmin)/100:xmax;
y = ymin:(ymax - ymin)/100:ymax;

map = zeros(length(y), length(x));

for i = 1:size(m, 2)
    norm = 1/sqrt((2*pi)^2 * det(V{i}(idx,idx)));
    prec = (V{i}(idx,idx))^-1;
    
    for j = 1:length(x)
        for k = 1:length(y)
            loc = [x(j) - m(idx(1), i); y(k) - m(idx(2), i)];
            map(k, j) = map(k, j) + W(i) * norm * exp( - 0.5 * loc' * prec * loc);
        end
    end
end

imagesc(x, y, map);
colormap hot
colorbar

map = map * (xmax - xmin)/100 * (ymax - ymin)/100;
end