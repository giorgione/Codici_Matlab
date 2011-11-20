%Energia Totole per il RMF

function Eloc=EnergiaLocale(i,j,X,Y,h,b,nn)

[m,n]=size(X);

%E2
%Indici dei vicini del pixel (i,j)
x=[-1 1 0  0];
y=[ 0 0 -1 1];
%E21=(X(i,j)*X(i,j+1)+X(i,j)*X(i,j-1)+X(i,j)*X(i-1,j)+X(i,j)*X(i-1,j));

E21=0;
for I=1:4
        i1=i+x(I);
        j1=j+y(I);
        if  (i1>0 && i1<=m) && (j1>0 && j1<=n)
            E21=E21+X(i,j)*X(i1,j1);
        end
end

E21=-b*E21;
E22=-nn*X(i,j)*Y(i,j);
E20=h*X(i,j);

Eloc=E20+E21+E22;



