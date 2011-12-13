function [Ic,Qc] = cropWindow( I, nBit, pos, wd, ht )
row = pos(2)-ht/2;  col = pos(1)-wd/2;
Ic = I(row:row+ht-1,col:col+wd-1,:);
if(nargout==2); Qc=bitshift(reshape(Ic,[],3),nBit-8); end;
end