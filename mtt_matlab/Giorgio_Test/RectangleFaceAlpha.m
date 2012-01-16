function RectangleFaceAlpha(Detections,FigHandle,Colore)
n=length(Detections);
figure(FigHandle)

for i=1:n
    % Get BBOX of Detection
    X=Detections(i).obs
    % X(1),X(2) --> Upper-Left corner
    % X(3) width
    % X(4) height
    if(isempty(X))
        break;
    end
    %
    w=X(3);
    h=X(4);
    
    %Upper-Left corner
    xo=X(1);
    yo=X(2);

    x1=xo+w;
    y1=yo;

    x2=x1;
    y2=y1+h;

    x3=xo;
    y3=y2;

    p=patch([xo x1 x2 x3],[yo y1 y2 y3],Colore);
    set(p,'FaceAlpha',0.2,'LineWidth',2);
end