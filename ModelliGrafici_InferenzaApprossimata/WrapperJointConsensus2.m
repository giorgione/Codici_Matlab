function P=WrapperJointConsensus2(Theta,M,N,Yij,a,b)
    Tau=Theta(1:M);
    u=Theta(M+1:end);
    P=JointConsensus(Tau,Yij,u,a,b);
end