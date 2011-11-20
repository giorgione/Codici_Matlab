function f = Lagrangiana_SVM(a,B)
[l,n]=size(a);
f=-(sum(a)-1/2*a.'*B*a);
