function E=EnergiaTotaleRMF_2(X,Y,h,b,n)
[m,n]=size(X);

%Indici dei vicini del pixel (i,j)
x=[ 1 0 ];
y=[ 0 1 ];

E1=h*(sum(sum(X)));

E2=0;
for i=1:m
    for j=1:n        
            for ngb=1:2
                i1=i+x(ngb);
                j1=j+y(ngb);
                if  (i1>0 && i1<=m) && (j1>0 && j1<=n)
                    E2=E2+X(i,j)*X(i1,j1);
                end
        end
    end
end
E2=-E2*b;

E3=0;
for i=1:m
    for j=1:n        
        E3=E3+X(i,j)*Y(i,j);
            
    end
end
E3=-E3*h;

E=E1+E2+E3;