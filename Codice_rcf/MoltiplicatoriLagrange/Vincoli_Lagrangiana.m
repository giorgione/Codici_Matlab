function [c, ceq] = Vincoli_Lagrangiana(a,y)
c = -a;
ceq = a.'*y;