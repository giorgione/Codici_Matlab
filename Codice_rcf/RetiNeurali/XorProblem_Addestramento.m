% XorProblem.m
% Risolvo il Problema dello Xor in R2 mediante una rete neurale in cui i le
% funzioni di Attivazione per l' HIDDEN LAYER sono di tipo signum (MLP)
%
% Idea:
%
% Ogni Livello della Rete Neurale viene mappato su di un Vettore.
clc;clear;close all
%Dati dello Xor Problem
X=[ 1 -1 -1  1 2 -2 -2  2 3 -3 -3  3 ;
    1 -1  1 -1 2 -2  2 -2 3 -3  3 -3 ];     
%Valori di Targhet
T=[ 0  0  1  1 0  0  1  1 0  0  1  1];
N=length(T);
n=N/2;
V=X;
V=[ones(1,N);V];


OutputLayer=@(x,a) Neurone([1;Neurone(x,[a(1);a(2);a(3)],'logis');Neurone(x,[a(4);a(5);a(6)],'logis')],[a(7);a(8);a(9)],'logis');

% Ho=[-3.9598 3.5405 1.7735;
%    -3.9598 -2.8277 2.772];
% Ho=Ho.';
% Oo=[-3.9531 4.1839 3.7222];
% ao=[Ho(:) ;Oo(:)];
ao=rand(9,1);

%Addestro la Rete mediante metodo del Gradiente discendente
options = optimset('Display','iter','TolFun',1e-5,'MaxIter',300,'TolX',1e-5);
Fun=@(a) ErroreGlobale(V,T,a);

for i=1:N
    Z(i)=OutputLayer(V(:,i),ao);
end
display('Output Iniziale')
disp(Z)


%Calcolo il minimo con fminsearch
[a1 ,Val,output]= fminsearch(Fun,ao,options);
for i=1:N
    Z(i)=OutputLayer(V(:,i),a1);
    if(Z(i)>.5)
        DisegnaClusters(V(2,i),V(3,i),1)
    else
        DisegnaClusters(V(2,i),V(3,i),-1)
    end
end
display('Algoritmo di Simulated Annealing: premere per continuare')
pause;
options=anneal();
options.InitTemp=1;
options.CoolSched=@(T)(.9*T);
options.Verbosity=2;
options.StopTemp=1e-16;
options.StopVal=1e-4;
options.MaxTries=500;
options.MaxConsRej=4500;
%options.Generator=@(x) (randn(length(x),1));
%Calcolo il Minimo mediante Simulated Annealing
[a2 f] = anneal(Fun,ao.',options);
figure;
for i=1:N
    Z(i)=OutputLayer(V(:,i),a2);
    if(Z(i)>.5)
        DisegnaClusters(V(2,i),V(3,i),1)
    else
        DisegnaClusters(V(2,i),V(3,i),-1)
    end
end
hold on;
% for i=1:N
%     Z(i)=OutputAddestrato(x(:,i));
%     %DisegnaClusters(x(2,i),x(3,i),Z(i))
% end

%Genero un Set di punti in [-1 1] x [-1 1]
% [x1 x2]=meshgrid(linspace(-2,2,20));
% [m,n]=size(x1);
% n=m*n;
% x=[ones(1,N); x1(:).'; x2(:).'];
% %Riclassifico e visualizzo i risultati
% figure;

