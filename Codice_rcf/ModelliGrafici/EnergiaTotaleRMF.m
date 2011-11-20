function E=EnergiaTotaleRMF(X,Y,h,b,nn)
[m,n]=size(X);

%Indici dei vicini del pixel (i,j)
x=[ -1 1  0 0 ];
y=[  0 0 -1 1 ];

E1=h*(sum(sum(X)));

E2=0;
for i=1:m
    for j=1:n 
        
            for ngb=1:4
                i1=i+x(ngb);
                j1=j+y(ngb);
                if  (i1>0 && i1<=m) && (j1>0 && j1<=n)
                    E2=E2+X(i,j)*X(i1,j1);
                end
            end
        
    end
end
E2=-b*E2;

E3=0;
for i=1:m
    for j=1:n        
        E3=E3+X(i,j)*Y(i,j);       
    end
end
E3=-nn*E3;

E=E1+E2+E3;


