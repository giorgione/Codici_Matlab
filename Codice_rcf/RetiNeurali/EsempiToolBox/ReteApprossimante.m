%Addestriamo una rete 

p = [0 1 2 3 4 5 6 7 8];
       t = [0 0.84 0.91 0.14 -0.77 -0.96 -0.28 0.66 0.99];
       plot(p,t,'o')
 
     %Here NEWFF is used to create a two layer feed forward network.
     %The network will have a single hidden layer of 10 neurons.
 
       net = newff(p,t,10);
       y1 = sim(net,p)
       plot(p,t,'o',p,y1,'x')
 
     %Here the network is trained for up to 50 epochs to a error goal of
     %0.01, and then resimulated.
 
       net.trainParam.epochs = 50;
       net.trainParam.goal = 0.01;
       net = train(net,p,t);
       y2 = sim(net,p)
       plot(p,t,'o',p,y1,'x',p,y2,'*')