function result=samstr(data,result);

D = data.X;
x = result.proj.P;

[noc,dim]=size(D);
  noc_x_1  = ones(noc, 1); 

  Mdist=[];
dim_x_1 = ones(dim,1);
  for i = 1:noc,
    xD = D(i,:); 
    Diff = D - xD(noc_x_1,:);
    Mdist(:,i) = sqrt((Diff.^2)*dim_x_1);
  end
  
  
  e = 0; tot = 0;
    for j = 2:noc, 
      d   = Mdist(1:(j - 1), j);
      tot = tot + sum(d);
      ind = find(d ~= 0);
      xd  = -x(1:(j - 1), :) + x(j * ones(j - 1, 1), :);
      ee  = d - sqrt(sum(xd'.^2))';
      e   = e + sum(ee(ind).^2 ./ d(ind));
    end
    e = e/tot; 
    
    result.proj.e = e;