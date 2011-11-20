function Pjoint=JointConsensus(Tau,Y,u,a,b)
    [M N]=size(Y);
    Pjoint=1;
    for j=1:M 
        Pjoint=Pjoint*gampdf(Tau(j),a,b);
        for i=1:N
               Pjoint=Pjoint*normpdf(Y(j,i),u,1/Tau(j));
        end
    end
end




