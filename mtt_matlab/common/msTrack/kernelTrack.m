%function [p, pos, Ic, sim] = kernelTrack( I, q, pos, kernel, nBit )
%
% Track the kernel inside the current Image I'
%
%
function [p, pos, Ic, sim] = kernelTrack( I, q, pos, kernel, nBit )
%Get image Dimensions
mRows=size(I,1); 
nCols=size(I,2);
%Get Kernel Dimensions
wd=kernel.wd;
wd2=wd/2;
ht=kernel.ht; 
ht2=ht/2;
%Get Kernel Samples in X and Y domain
xs=kernel.xs; %Vector
ys=kernel.ys; %Vector

for iter=1:1000
      %Get Previous Position
      posPrev = pos;

      % check if pos in bounds
      % calculate the rectangle (wd x ht) centered in pos = [uc vc w h]
      rct = [pos(1)-wd/2 pos(2)-ht/2 wd, ht ];
      %check if the rectangle is into the image 
      if( rct(1)<1 || rct(2)<1 || (rct(1)+wd)>nCols || (rct(2)+ht)>mRows )
        %in this case we have an invalid position so we have the solution
        pos=posPrev; 
        p=[];
        Ic=[]; 
        sim=eps; 
        return;
      end

      % crop window / compute histogram on the current frame
      [Ic,Qc] = cropWindow( I, nBit, int16(pos), wd, ht );
     
      %build histogram on the Window: p
      p = buildHist( Qc, kernel, nBit );
      if( iter==20 );
          break;
      end;

      % compute meanshift step -- Mean Shift Weight
      w = ktComputeW_c( Qc, q, p, nBit );

      posDel = [sum(xs.*w)*wd2, sum(ys.*w)*ht2] / (sum(w)+eps);
      posDel = round(posDel+.1);
      if( all(posDel==0) ) 
          break; 
      end

      pos = pos + posDel;
end
%Find positive bins in the Histogram p
locs=p>0;

%Calculate the  Similiarity Measure
sim=sum( sqrt(q(locs).*p(locs)) );

end