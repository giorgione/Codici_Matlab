function showKLT( dirname, ext, kltfile )
[x, y, val] = klt_read_featuretable(kltfile);

files = dir([dirname '/*.' ext]);
for i = 1:length(files)
    im = imread([dirname '/' files(i).name]);
    imshow(im);
    hold on;
    scatter(x(:, i), y(:, i), 'b.');
    if i > 1
        tidx = find(val(:, i) == 0);
        quiver(x(tidx, i-1), y(tidx, i-1), x(tidx, i) - x(tidx, i - 1), y(tidx, i) - y(tidx, i - 1), .4, 'r');
    end
    hold off;
    pause(0.5);
end

end
