function content = parseSignAnnotations(fileName)
%function content = parseSignAnnotations(fileName)
%
% Given a file name, the function parses the corresponding text file and
% returns a struct with sign info.
%
% e.g.
% content = parseSignAnnotations('annotations.txt')
%
%/larsson@isy.liu.se 2011

fid = fopen(fileName);
% h = figure(1);

imageNr = 0;
oo=0;
while ~feof(fid)
	
	line = fgetl(fid);
	nameEnd = find(line == ':');
	if nameEnd == []
		disp('ERROR: file incorrect!');
		break;
	end
	
	imName = line(1:nameEnd-1);
			
	imageNr = imageNr + 1;
	
	content(imageNr).name              = imName;
	content(imageNr).signs.signTypes   = [];
	content(imageNr).signs.signBB      = [];
	content(imageNr).signs.signC       = [];
	content(imageNr).signs.signStatus  = [];
	content(imageNr).signs.signSize    = [];
	content(imageNr).signs.aspectRatio = [];
	
			
	if length(line) == nameEnd
		%No sign in this image;
	else
		nrSigns = 0;
		signEnd = find(line == ';');
		while length(signEnd > 0)
		%misc. signs found in image
		if isequal(line(nameEnd+1:signEnd(1)-1),'MISC_SIGNS')
			nrSigns = nrSigns + 1;
			line = line(signEnd(1)+1:end);
			signEnd = find(line == ';');
			
			content(imageNr).signs(nrSigns).signTypes   = 'MISC_SIGNS';
			content(imageNr).signs(nrSigns).signStatus  = 'N/A';
			content(imageNr).signs(nrSigns).signBB      = [-1 -1 -1 -1];			
			content(imageNr).signs(nrSigns).signC       = [-1 -1];			
			content(imageNr).signs(nrSigns).signSize    = 0;			
			content(imageNr).signs(nrSigns).aspectRatio = 0;
			
		else						
			
				nrSigns = nrSigns + 1;
				commas = find(line == ',');
				
				%Extract sign information
				visibility = line(1:commas(1)-1);
 				%remove image name if present
 				removeIndex = strfind(visibility,':');
				if ~isempty(removeIndex)
					visibility = visibility(removeIndex+1:end);
				end
				
				lrx = str2double(line(commas(1)+2:commas(2)-1));
				lry = str2double(line(commas(2)+2:commas(3)-1));
				ulx = str2double(line(commas(3)+2:commas(4)-1));
				uly = str2double(line(commas(4)+2:commas(5)-1));
				signType = line(commas(5)+2:commas(6)-1);
				signName = line(commas(6)+2:signEnd(1)-1);
							
				line = line(signEnd(1)+1:end);
				signEnd = find(line == ';');
												
				content(imageNr).signs(nrSigns).signTypes   = signName;
				content(imageNr).signs(nrSigns).signStatus  = visibility;
				content(imageNr).signs(nrSigns).signBB      = [ulx uly lrx lry];
				content(imageNr).signs(nrSigns).signC       = [(ulx+lrx)/2 (uly+lry)/2];
				content(imageNr).signs(nrSigns).signSize    = min((ulx-lrx),(uly-lry))^2;
				content(imageNr).signs(nrSigns).aspectRatio = max((ulx-lrx)/((uly-lry)+eps),(uly-lry)/((ulx-lrx)+eps)); 
				
				
			end
		end
	end	
end
fclose(fid);

