function nrSigns = displaySigns(signInfo,options)
% function nrSigns = displaySigns(signInfo,options)
% Displays the signs given in a signInfo-struct.
%
% options = 0 : Only displays bounding boxes around signs
% options = 1 : Displays bounding box and sign type
% options = 2 : Displays bounding box, sign type and additional info
%
% /larsson@isy.liu.se 2011



nrSigns = size(signInfo,2);
if nrSigns>0
	if(isempty(signInfo(1).signTypes))
		nrSigns = 0;
		return;
	end
end

if nargin<2
	options=1;
end


%Plot parameters
lw = 2;   %line width
lc = 'r'; %line color
textOffset = 25;


for i=1:nrSigns
	info = signInfo(i);
	
	if strcmpi(info.signTypes,'MISC_SIGNS')
		if options>0
			%Write text in upper left corner
			h=text(1,20,escape_underscores(info.signTypes),'color',lc);
			%set(h,'BackgroundColor',[0 0 0]);
		end
	else
		%plot bounding box
		bb = info.signBB;
		x = [bb(1) bb(1) bb(3) bb(3) bb(1)];
		y = [bb(2) bb(4) bb(4) bb(2) bb(2)];
		hold on;
		plot(x,y,lc,'lineWidth',lw);
		hold off;
		
		if options>0
			%Add text
			x = max(bb(1)-textOffset,1);
			y = max(bb(2)-textOffset,1);
			if options<2
				h=text(x,y,escape_underscores(info.signTypes),'color',lc);
			else
				h=text(x,y,escape_underscores([info.signTypes,', ',info.signStatus]),'color',lc);
			end
			%set(h,'BackgroundColor',[0 0 0]);
		end
	end
end







%function s=escape_underscores(s0)
%
% Convert a string to a new string where
% underscores are converted to '\_' for
% use in Matlab display strings.
%
function s=escape_underscores(s0)

upos=find(s0=='_');
Nu=numel(upos);

if ~isempty(upos)
	
	s=[s0(1:(upos(1)-1)) '\_'];
	
	for k=2:Nu,
		s=[s s0(upos(k-1)+1:upos(k)-1) '\_'];
	end
	s=[s s0(upos(end)+1:end)];
else
	s=s0;
end