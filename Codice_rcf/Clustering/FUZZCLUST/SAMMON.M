function result = Sammon(proj,data,result,param)
%function P = sammon(D, P, varargin)

%SAMMON Computes Sammon's mapping of a data set.


%  Input and output arguments ([]'s are optional):
%   D        (matrix) size dlen x dim, data to be projected
%            (struct) data or map struct            
%   P        (scalar) output dimension
%            (matrix) size dlen x odim, initial projection matrix
%   [value]  (scalar) all different modes (the next argument) require 
%                     a value, default = 100
%   [mode]   (string) 'steps' or 'errlimit' or 'errchange' or 'seconds',
%                     see below, default is 'steps'
%   [alpha]  (scalar) iteration step size, default = 0.2
%   [Dist]   (matrix) pairwise distance matrix, size dlen x dlen.
%
%   P        (matrix) size dlen x odim, the projections
%
% The output dimension must be 2 or higher but (naturally) lower 
% than data set dimension.
%
% The mode argument determines the end condition for iteration. If 
% the mode argument is used, also the value argument has to be 
% specified. Different mode possibilities are:
% 'steps'      the iteration is terminated when it is run <value> 
% 'errlimit'   steps, the iteration is terminated when projection error 
%              is lower than <value>,
% 'errchange'  the iteration is terminated when change between 
%              projection error on two successive iteration rounds
%	       is less than <value> percent of total error, and
% 'seconds'    the iteration is terminated after <value> seconds 
%              of iteration.
%

% Reference: Sammon, J.W. Jr., "A nonlinear mapping for data
%   structure analysis", IEEE Transactions on Computers, vol. C-18,
%   no. 5, 1969, pp. 401-409.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check arguments
%Adapt the data
P = proj.P;             %P projected data or dimension
D = data.X;
%Mdist = result.data.d;      %d distances
maxstep = param.max;
alpha = param.alpha;
m = param.m;

% compute data dimensions
orig_si = size(D); 
dim = orig_si(2); 
noc = orig_si(1);

% output dimension / initial projection matrix
if prod(size(P))==1, 
  odim = P; 
  P = rand(noc,odim)-0.5; 
else 
  si = size(P);
  odim = si(end);
  if prod(si) ~= noc*odim, 
    error('Initial projection matrix size does not match data size');
  end
end
if odim > dim | odim < 2, 
  error('Output dimension must be within [2, dimension of data]');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialization

% these are used quite frequently
noc_x_1  = ones(noc, 1); 
odim_x_1 = ones(odim,1); 

%compute mutual distances between vectors
Mdist = [];
if isempty(Mdist) | all(isnan(Mdist(:))),  
  fprintf(2, 'computing mutual distances\r');
  dim_x_1 = ones(dim,1);
  for i = 1:noc,
    x = D(i,:); 
    Diff = D - x(noc_x_1,:);
    Mdist(:,i) = sqrt((Diff.^2)*dim_x_1);
  end
end
[np,nc]=size(Mdist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% action
% sammon iteration   

x  = P ;
xu = zeros(noc, odim);
xd = zeros(noc, odim);
dq = zeros(noc, 1);
dr = zeros(noc, 1);

for i=1:maxstep
    
   % If you want to see the Sammon's projection plotted (in 2-D and 3-D case),
  % execute the code below; it is not in use by default to speed up 
  % computation.  
  if 1, 
    clf
    if odim == 1,     plot(x(:,1), noc_x_1, 'o');
    elseif odim == 2, plot(x(:,1), x(:,2), 'o');
    else              plot3(x(:,1), x(:,2), x(:,3), 'o')
    end
    drawnow
  end
  
  
   
    
    
  for j = 1:noc,
    xd      = -x + x(j*noc_x_1,:);
    xd2     = xd.^2;
    dpj     = sqrt(sum(xd2'))';
    dq      = Mdist(:,j) - dpj;
    dr      = Mdist(:,j) .* dpj;
    ind     = find(dr ~= 0);
    term    = dq(ind) ./ dr(ind);
    e1      = sum(xd(ind,:) .* term(:,odim_x_1));
    term2   = ((1.0 + dq(ind) ./ dpj(ind)) ./ dpj(ind)) ./ dr(ind);
    e2      = sum(term) - sum(xd2(ind,:) .* term2(:,odim_x_1));
    xu(j,:) = x(j,:) + alpha * e1 ./ abs(e2);
  end

  % move the center of mass to the center 

  c = sum(xu) / noc;
  x = xu - c(noc_x_1, :);

  
  fprintf(2, '\r%d iterations', i);
  
  
  
  fprintf(2, '        ');
  
  
end

fprintf(2, '\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clean up

% reshape
orig_si(end) = odim; 
P = reshape(x, orig_si);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%computing cluster centers
fm = result.data.f.^m;
sumf = sum(fm);
vs = (fm'*x)./(sumf'*ones(1,2));

%%%%%%%%%%%%%%%%%%%%%%
result.proj.vp=vs;
result.proj.P=x;