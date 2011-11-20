syms s1 u1 P1 s2 u2 P2 x R;


ps1=GaussianaMulti(s1,u1,P1);
ps2=GaussianaMulti(s1,u1,P2);

px_s1s2=GaussianaMulti(x,s1+s2,R);

%Il modello grafico  mi dice che la Prob. Congiunta p(s1,s2,x)
PX=ps1*ps2*px_s1s2;

PX=subs(PX,{u1,u2,P1,P2,R},{3,5,0.5,0.5,.3})
Px=int(int(PX,s1),s2);

Ps1s2_x=PX/Px;
Ps1s2_x=subs(Ps1s2_x,x,9)






