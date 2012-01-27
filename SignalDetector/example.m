%Script that reads a content files and displays annotation info
%/larsson@isy.liu.se 2011

%Set paths to the images and to the annotation file
base = './';
imagePath = fullfile(base,'Set2Part0');
annotationFn = fullfile(imagePath,'annotations.txt');

 %Get content struct from annotation file
content = parseSignAnnotations(annotationFn);      
content1=EstraiSegnali(content,'70_SIGN')
%For all images with annotation, display content
N = length(content);

for i = 1:N
	fn = fullfile(imagePath,content(i).name);
	try
		img = imread(fn);
		imagesc(img);axis image;
		title(content1(i).name);
	catch
		display(['Could not read image: ',fn]);		
	end
	
	try
		%Display content (only bounding boxes)
		displaySigns(content1(i).signs,2 );
	catch
		display('could not display sign info, something is wrong');
    end	
    
	pause
    %content1(i).signs.name
end
