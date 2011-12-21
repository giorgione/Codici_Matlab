%% CAMERA initial paramters:  Variable 8x1 parameters
% focal length, camera height(meter), x camera center(pixel), horizon
% (pixel), initial angle, initial speed, camera's x location, camera's z location
%Z.cam = [1200 1.0 imsize(1)/2 ihorizon 0 11 0 0]';

function DrawSample(sample,myFig,Colore)


figure(myFig)
%Assi Mondo
Aw=100*eye(3);
%Assi Mondo
plot3(0,0,0,'or','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5);
coloreM=[1 0 0];
colorX=[1 0 0];
colorY=[0 1 0];
colorZ=[0 0 1];
%Sistema di Riferimento Mondo centrato nell'origine
Vettore3D(Aw(1,1),Aw(2,1),Aw(3,1),colorX,myFig) %X in red
Vettore3D(Aw(1,2),Aw(2,2),Aw(3,2),colorY,myFig) %Z in green
Vettore3D(Aw(1,3),Aw(2,3),Aw(3,3),colorZ,myFig) %Y in blue

%Disegno la Camera 
Cam=sample.cam;
[m,n]=size(Cam);
for i=1:n
    
    K=@(f,ku,kc,uc,vc)[f*ku 0 uc
                       0 f*kv vc
                       0 0  1];



    Rc2w =@(tx) [cos(tx)  sin(tx)  0; 
                 -sin(tx)  cos(tx)  0;
                    0          0   1];

    Rw2c=@(tx)inv(Rc2w(tx));

    %World to Camera
    C2W=@(tx,Zc,T) Rc2w(tx)*(Zc-T);

    %Camera to World
    W2C=@(tx,Zw,T) Rw2c(tx)*Zw+T;

    ProJ=@(Zw,f,ku,kv,uc,vc,tx,T) K(f,ku,kv,uc,vc )*W2C(tx,Zw,T) ;
    %Fisso i PARAMETRI INTRINSECI
    fo=Cam(1,i);
    Uo=Cam(3,i);
    Vo=Cam(3,i);
    
    ku=10^-6; kv=10^-6;
    w=2*Cam(3,i)*ku;
    h=2*Cam(3,i)*kv; 
    %PARAMETRI ESTRINSECI
    tx=Cam(2);

    % Posizione del Centro Ottico C0
    Xt=Cam(7,i);
    Zt=Cam(8,i); %profondita
    Ht=0;  %altezza della camera
    T=[Xt Zt Ht].';
    CO=T;
    
    plot3(T(1),T(2),T(3),'o','MarkerFaceColor',Colore,'MarkerEdgeColor',Colore,'MarkerSize',5);
    


    %Fissato un angolo di Rotazione vado a vedere la relazione che intercorre
    %trai 2 sistemi di riferimento
    %Verifico che RW2C e' una ROTAZIONE sull' asse Z
    Ac=W2C(tx,Aw,repmat([0;0;0],1,3));
    %Sistema di Riferimento Mondo centrato nell'origine
    %Vettore3D(Ac(1,1),Ac(2,1),Ac(3,1),colorX,1) %X in red
    %Vettore3D(Ac(1,2),Ac(2,2),Ac(3,2),colorY,1) %Y in green
    %Vettore3D(Ac(1,3),Ac(2,3),Ac(3,3),colorZ,1) %Z in blue

    %Disegno i Vettori applicati in T
    %La camera e allineata al piano di terra e l'asse focale e l'asse Y
    AsseU=w*Ac(:,1)./norm(Ac(:,1));
    AsseV=h*Ac(:,3)./norm(Ac(:,3));
    DisegnaPiano(AsseU,AsseV,T,myFig,Colore);
    Ac=Ac/10;
    Vettore3D_Applicato(T,Ac(1,1),Ac(2,1),Ac(3,1),colorX,myFig)
    Vettore3D_Applicato(T,Ac(1,2),Ac(2,2),Ac(3,2),colorY,myFig)
    Vettore3D_Applicato(T,Ac(1,3),Ac(2,3),Ac(3,3),Colore,myFig)
    drawnow
end
Person=sample.per;
[m,n]=size(Person);
for i=1:n
    if(m==0)
        break;
    end
    %Disegno la persona
    %P=C2W(tx,Person(1:3,i),T);
    plot3(Person(1,i),Person(2,i),Person(3,i),'o','MarkerFaceColor',Colore,'MarkerEdgeColor',Colore,'MarkerSize',5);
end

Car=sample.car;
[m,n]=size(Car);
for i=1:n
    if(m==0)
        break;
    end
    %Disegno la MACCHINA
    %C=C2W(tx,Car(1:3,i),T);
    plot3(Car(1,i),Car(2,i),Car(3,i),'o','MarkerFaceColor',Colore,'MarkerEdgeColor',Colore,'MarkerSize',5);

end
