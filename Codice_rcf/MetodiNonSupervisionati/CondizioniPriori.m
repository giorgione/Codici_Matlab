function [c ceq]=CondizioniPriori(x)
% Nonlinear inequality constraints
c = [-x(5);     
     -x(6)];
% Nonlinear equality constraints
ceq = [x(5)+x(6)-1];