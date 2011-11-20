clc;clear;close all


%Test sulle Reti Neurali (Problema dello XOR) con il TOOLBOX di Matlab
%funzioni di Attivazione TANSIG
Pattern=[-1 1  1 -1; 
         -1 1 -1  1];
     
Target=[-1 -1 1 1];

net=newff(Pattern,Target,2,{},'traingdm');
net.trainParam.show=50;
net.trainParam.lr=0.05;
net.trainParam.mc = 0.9;
net.trainParam.epochs=500;
net.trainParam.goal=1e-5;

[net,tr]=train(net,Pattern,Target); % tr (training record) contiene informazioni su come
% procede l'addestramento
Out=sim(net,Pattern)
figure;
for i=1:4
    DisegnaClusters(Pattern(1,i),Pattern(2,i),sign(Out(i)))
end