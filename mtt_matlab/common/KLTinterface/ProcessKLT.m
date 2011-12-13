function ProcessKLT(dirname, ext, nfeat)
files = dir([dirname '/*.' ext]);

outformat = 'temp/img%06d.pgm';
for i = 1:length(files)
    outfile = sprintf(outformat, i - 1);
    if(exist(outfile))
        continue;
    end
    try
        im = imread([dirname '/' files(i).name]);
    catch
        files(i+1:end) = [];
        break;
    end
    imwrite(im, outfile, 'pgm');
end

system(['../../../KLT ' num2str(length(files)) ' ' num2str(nfeat)]);

for i = 1:length(files)
    outfile = sprintf(outformat, i);
    delete(outfile);
end

movefile('features.txt', [dirname '/features.txt']);

end
