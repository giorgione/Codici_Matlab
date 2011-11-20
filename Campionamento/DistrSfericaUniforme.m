function p=DistrSfericaUniforme(x,r,d)

%Volume della Sfera
Vol=pi^(d/2) *( r^2)/ gamma(d/2);

if x'*x <= r^2
    p=1/Volume;
else
    p=0;