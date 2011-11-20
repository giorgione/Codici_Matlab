function [H]=hessian(f,vars)
%% function [H]=hessian(f,vars)
% Computes the Hessian of a given function symbolically
% Usage Example:
%     >> syms x1 x2 x3;
%     >> f=x1^4+x2^4+x3^2;
%     >> hessian(f,[x1 x2 x3])
%  
%     ans =
%  
%        [ 12*x1^2,       0,       0]
%        [       0, 12*x2^2,       0]
%        [       0,       0,       2]

numVars=length(vars);
for(row=1:numVars)
    for(col=1:numVars)
        H(row,col)=diff(diff(f,vars(row)),vars(col));
    end
end

