% [X,F] = Obfgs2(f,paramf,x0[,tol,maxit,qN,H0]) Optimisation by BFGS
%
% Reference: p. 194ff in:
%   Nocedal and Wright: "Numerical Optimization", Springer, 1999.
%
% Example usage & In/Out arguments: see Osteepdesc.m.
%
% In:
%   H0: n x n matrix, initial inverse Hessian approximation. Default: eye(n,n).
%   qN: quasi-Newton method, one of 'BFGS', 'DFP', 'SR1'. Default: 'BFGS'.
%
% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2006 by Miguel A. Carreira-Perpinan

function [X,F] = Obfgs2(f,paramf,x0,tol,maxit,qN,H0)

n = length(x0);		% Dimension of the problem.

% ---------- Argument defaults ----------
if ~exist('tol','var') || isempty(tol) tol = 1e-5; end;
if ~exist('maxit','var') || isempty(maxit) maxit = 5000; end;
if ~exist('qN','var') || isempty(qN) qN = 'BFGS'; end;
if ~exist('H0','var') || isempty(H0) H0 = eye(n,n); end;
% ---------- End of "argument defaults" ----------

x = x0;
[ff,g] = f(x',paramf{:});
g = g';
X = x';
F = ff;
H = H0;
k = 0;
ng = 0;

oldx = zeros(size(x));         % add by cqin

%fprintf(1,'Cycle %d: ff=%f norm(g)=%f\n',k,ff,norm(g));

cnt = 0; max_cnt = 5;

%while (norm(g) > tol) & (k < maxit)  % original version
while (norm(g) > tol) & (norm(x-oldx) > 1e-4) & (k < maxit) 
  oldx = x;
  oldg = g;
  
  p = -H*g;
  if p'*g>=0
    p = -g;
    fprintf('Change to negative gradient!\n');
    ng = ng+1;
  end
  
  % Linesearch to satisfy Wolfe condition
  alpha = stronglinesearch(f,paramf,x,p,ff,g',1);
  
  s = alpha*p;
  x = x + s;				% Get the new x.
  
  [ff,g] = f(x',paramf{:});		% Get the new gradient.
  g = g';
  y = g - oldg;
  rho = 1/(y'*s);
  if strcmp(upper(qN),'BFGS')
    if y'*s>0
      H = (eye(n,n)-rho*s*y')*H*(eye(n,n)-rho*y*s')+rho*(s*s');	% BFGS update.
    end
  else
    if strcmp(upper(qN),'DFP')
      H = H-(H*y)*(y'*H)/(y'*H*y)+rho*(s*s');		% DFP update.
    else
      H = H+(s-H*y)*(s-H*y)'/((s-H*y)'*y);		% SR1 update.
    end
  end
  
  %if norm(g)==norm(oldg)  % original
  if (abs(norm(g)-norm(oldg))<1e-5)
    cnt = cnt + 1;
    if (cnt>max_cnt) break; end;
  else
    cnt = 0;    % reset
  end
  
  X = [X;x'];
  F = [F;ff];
  k = k + 1;
  
  %fprintf(1,'Cycle %d: ff=%f norm(g)=%f norm(dx)=%f alpha=%f r=%f\n',k,ff,norm(g),norm(x-oldx),alpha,y'*s);
  
end

%fprintf('Change to negative gradient for totoally %d times!\n',ng)

if (k==maxit)
  fprintf(1,'Max iterations reached!\n');
end
