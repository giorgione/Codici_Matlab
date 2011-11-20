function Pn=SommaRettangoli(X,x,h,Vn)
n=length(x);

Pn=0;
for i=1:n
    Pn=Pn+rettangolo((X-x(i))/h,h);
end

Pn=Pn/Vn;