function kernel = buildKernel( wd, ht )
wd = round(wd/2)*2;  xs = linspace(-1,1,wd);
ht = round(ht/2)*2;  ys = linspace(-1,1,ht);
[ys,xs] = ndgrid(ys,xs); xs=xs(:); ys=ys(:);
xMag = ys.*ys + xs.*xs;  xMag(xMag>1) = 1;
K = 2/pi * (1-xMag);  sumK=sum(K);
kernel = struct( 'K',K, 'sumK',sumK, 'xs',xs, 'ys',ys, 'wd',wd, 'ht',ht );
end
