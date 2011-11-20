% function DisegnaClusters3D(X,Y,F)
% Disegna i punto (X,Y) con il colore specificati dal vettore F
%
% parametri in:
%
% X: ascissa del punto
%
% Y: ordinata del punto
%
% F: colore del punto i-esimo
%                               1 - blue
%                               2 - verde
%                               3 - rosso
%                               4 - magenta

function DisegnaClusters3D(X,Y,Z,F)
size=length(F);
for i=1:size
    switch F(i)
        case 1
            color='b';
        case 2
            color='g';
        case 3
            color='r';
        case 4
            color='m';
    end
            
plot3(X(i),Y(i),Z(i),['o' color],'MarkerSize',3,'MarkerFaceColor',color); hold on;

end
