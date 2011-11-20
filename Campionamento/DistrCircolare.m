function P=DistrCircolare(X,Xc,R)
if((X-Xc)'*(X-Xc) <= R^2)
    P=1/((R^2)*pi);
else
    P=0.0;
end