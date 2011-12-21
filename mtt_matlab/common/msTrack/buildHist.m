function p = buildHist( Qc, kernel, nBit )
p = ktHistcRgb_c( Qc, kernel.K, nBit ) / kernel.sumK;
if(0);p=gaussSmooth(p,.5,'same',2); p=p*(1/sum(p(:))); end;
end
