function hist = dbgComputeColHist(img, colStep)

imHeight = size(img, 1);
imWidth = size(img, 2);

colStep = 0-colStep;
colSize = bitshift(256, colStep);

hist = zeros(colSize, colSize, colSize);

for i = 1:imHeight
    for j = 1:imWidth
        idx1 = bitshift(img(i, j, 1), colStep) + 1;
        idx2 = bitshift(img(i, j, 2), colStep) + 1;
        idx3 = bitshift(img(i, j, 3), colStep) + 1;
        
        hist(idx1, idx2, idx3) = hist(idx1, idx2, idx3) + 1/(imWidth*imHeight);
    end
end

end