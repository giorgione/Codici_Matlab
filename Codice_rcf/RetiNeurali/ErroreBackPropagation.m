function E=ErroreBackPropagation(Pattern,Target,a)
[Dim NPattern]=size(Pattern);
E=0;
for i=1:NPattern
    Yk=Neurone( [Neurone(Pattern(:,i),[a(1);a(2)],'logis');Neurone(Pattern(:,i),[a(3);a(4)],'logis')],[a(5);a(6)],'logis');
    E=E+1/2*(Target(i)-Yk)^2;
end
