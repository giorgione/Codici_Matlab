function [p, pos, Ic, sim] = kernelTrack( I, q, pos, kernel, nBit )

mRows=size(I,1); nCols=size(I,2);
wd=kernel.wd; wd2=wd/2;
ht=kernel.ht; ht2=ht/2;
xs=kernel.xs; ys=kernel.ys;

for iter=1:1000
  posPrev = pos;

  % check if pos in bounds
  rct = [pos(1)-wd/2 pos(2)-ht/2 wd, ht ];
  if( rct(1)<1 || rct(2)<1 || (rct(1)+wd)>nCols || (rct(2)+ht)>mRows )
    pos=posPrev; p=[]; Ic=[]; sim=eps; return;
  end

  % crop window / compute histogram
  [Ic,Qc] = cropWindow( I, nBit, int16(pos), wd, ht );
  p = buildHist( Qc, kernel, nBit );
  if( iter==20 ); break; end;

  % compute meanshift step
  w = ktComputeW_c( Qc, q, p, nBit );
  posDel = [sum(xs.*w)*wd2, sum(ys.*w)*ht2] / (sum(w)+eps);
  posDel = round(posDel+.1);
  if(all(posDel==0)); break; end;
  pos = pos + posDel;
end

locs=p>0; sim=sum( sqrt(q(locs).*p(locs)) );
end