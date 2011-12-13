function d=bhattacharyya(X1,X2)

bc = sum(sqrt(X1.*X2));
d = -log(bc);

end
