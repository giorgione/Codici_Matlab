% function S=MatriceScatter(x,Media)
% Calcolo la Matrice di Scatter associata alle osservazioni x e al vettore
% Media
function S=MatriceScatter(x,Media)
    [m,n]=size(x)
    S=zeros(m,m);
    for i=1:n
        S=S+(x(:,i)-Media)*(x(:,i)-Media).';
    end