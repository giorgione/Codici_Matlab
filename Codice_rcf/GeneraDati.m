function Data=GeneraDati(dim,N,type)


switch lower(type)
          case {'s'}
            %genera i dati SEPARATI
            C1=10*rand(dim,N);
            C2=200+10*rand(dim,N);
            Data=[C1 C2];
          case 'r'
            %genera i dati NON SEPARATI
            Data=rand(dim,2*N)
end