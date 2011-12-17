% Draw 3D predicion
function Draw3DTrackingResult(Myfig,X,type)
if type==1
    %CAR
    colore='r'
else
    %PERSON
    colore='b'
end

figure(Myfig);hold on

plot3(X(1),X(2),X(3),'oy','MarkerFaceColor',colore,'MarkerEdgeColor',colore,'MarkerSize',5);
