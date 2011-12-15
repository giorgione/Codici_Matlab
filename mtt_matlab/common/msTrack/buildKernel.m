%function kernel = buildKernel( wd, ht )
%Generate the RECTANGULAR KERNEL of wd x ht points in [-1 1]x[-1 1] for the MEAN SHIFT TRAKER
% 
%               2
% K(X)  = ------------------
%                          2
%         pi* ( 1 - || x || )
%
function kernel = buildKernel( wd, ht )
    %Make exactly divisble by 2 the number of points
    wd = round(wd/2)*2; 
    ht = round(ht/2)*2;  
    %generate  wd x ht points in [-1 1]x[-1 1] 
    xs = linspace(-1,1,wd);
    ys = linspace(-1,1,ht);
    [ys,xs] = ndgrid(ys,xs); 
    xs=xs(:); ys=ys(:);
    
    %Module
    xMag = ys.*ys + xs.*xs;  
    
    xMag(xMag>1) = 1;
    K = 2/pi * (1-xMag);
    %Sum
    sumK=sum(K);
    kernel = struct( 'K',K, 'sumK',sumK, 'xs',xs, 'ys',ys, 'wd',wd, 'ht',ht );
end
