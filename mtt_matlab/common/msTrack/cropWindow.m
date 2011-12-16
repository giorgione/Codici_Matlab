%function [Ic,Qc] = cropWindow( I, nBit, pos, wd, ht )
% Crop a window Ic(wd x ht) centered in pos from the Image I and save the
% result of bitshift with nBit-8 position in Qc
%
% PARAMETERS INPUT:
%
%   - I:Input Image 
%   - nBit: Number of Bit for shifting 
%   - pos: Window Center in the Image
%   - wd, ht: Window dimensions
%
% PARAMETERS OUT_:
%  - Ic: the extracted Window from I
%  - Qc: Ic bitshifted with nBit-8
function [Ic,Qc] = cropWindow( I, nBit, pos, wd, ht )
    %Extract the Window
    row = pos(2)-ht/2;  
    col = pos(1)-wd/2;
    Ic = I(row:row+ht-1,col:col+wd-1,:);
    %Save in Qc the result of the Bit Shift Operation on Ic
    if(nargout==2); Qc=bitshift(reshape(Ic,[],3),nBit-8); end;
end