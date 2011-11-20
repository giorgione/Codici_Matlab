function Pn=SommaRettangoli2(X,x,h)
n=length(x);

Pn=0;
for i=1:n
    Pn=Pn+rettangolo((X-x(i))/h,1);
end

Pn=Pn/h;