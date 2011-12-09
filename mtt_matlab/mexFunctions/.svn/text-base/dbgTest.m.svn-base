% test image similarity
clear

im1 = imread('imgs/ex1.jpg');
hist1 = mexComputeAppearanceLkhood(im1, [], 4, 1, 3);

im2 = imread('imgs/ex2.jpg');
hist2 = mexComputeAppearanceLkhood(im2, [], 4, 1, 3);

im3 = imread('imgs/ex3.jpg');
hist3 = mexComputeAppearanceLkhood(im3, [], 4, 1, 3);

lkhood1 = zeros(1, 10);
for test = 1:10
    J1 = imnoise(im1,'salt & pepper', 0.01 * test);
    [temp, lkhood1(test)]=mexComputeAppearanceLkhood(J1, hist1, 4, 1, 3);
end


lkhood2 = zeros(1, 10);
for test = 1:10
    J2 = imnoise(im2,'salt & pepper', 0.01 * test);
    [temp, lkhood2(test)]=mexComputeAppearanceLkhood(J2, hist2, 4, 1, 3);
end


lkhood3 = zeros(1, 10);
for test = 1:10
    J3 = imnoise(im3,'salt & pepper', 0.01 * test);
    [temp, lkhood3(test)] = mexComputeAppearanceLkhood(J3, hist3, 4, 1, 3);
end


lkhood12 = zeros(1, 10);
for test = 1:10
    J2 = imnoise(im2,'salt & pepper', 0.01 * test);
    [temp, lkhood12(test)]=mexComputeAppearanceLkhood(J2, hist1, 4, 1, 3);
end


lkhood23 = zeros(1, 10);
for test = 1:10
    J3 = imnoise(im3,'salt & pepper', 0.01 * test);
    [temp, lkhood23(test)]=mexComputeAppearanceLkhood(J3, hist2, 4, 1, 3);
end

lkhood13 = zeros(1, 10);
for test = 1:10
    J3 = imnoise(im3,'salt & pepper', 0.01 * test);
    [temp, lkhood13(test)]=mexComputeAppearanceLkhood(J3, hist1, 4, 1, 3);
end


lkhood = [];
imf = imread('imgs/exFull.jpg');
for i = 1:(size(imf, 1) - size(im1, 1))
    for j = 1:(size(imf, 2) - size(im1, 2))
        tempim = imf(i:(i+size(im1, 1)-1), j:(j+size(im1, 2)-1), :);
        [temp, lkhood(i, j)]=mexComputeAppearanceLkhood(tempim, hist1, 4, 1, 3);
    end
end
