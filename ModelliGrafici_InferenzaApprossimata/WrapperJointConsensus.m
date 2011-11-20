function P=WrapperJointConsensus(Theta,M,N,u,a,b)
    Tau=Theta(1:M);
    y=Theta(M+1:end);
    Yij=reshape(y,M,N);
    P=JointConsensus(Tau,Yij,u,a,b);
end

