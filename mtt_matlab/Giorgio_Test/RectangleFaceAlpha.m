function RectangleFaceAlpha(X,FigHandle)

w=X(3);h=X(4)
%Upper-Left corner
xo=X(1);
yo=X(2);

x1=xo+w;
y1=yo;

x2=x1;
y2=y1+h;

x3=xo;
y3=y2;
figure(FigHandle)
p=patch([xo x1 x2 x3],[yo y1 y2 y3],'r');
set(p,'FaceAlpha',0.2); 