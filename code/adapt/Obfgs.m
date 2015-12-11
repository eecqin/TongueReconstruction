% [X,F] = Obfgs(f,paramf,x0[,tol,maxit,qN,H0]) Optimisation by BFGS
%
% Reference: p. 194ff in:
%   Nocedal and Wright: "Numerical Optimization", Springer, 1999.
%
% Example usage & In/Out arguments: see Osteepdesc.m.
%
% In:
%   H0: n x n matrix, initial inverse Hessian approximation. Default: eye(n).
%   qN: quasi-Newton method, one of 'BFGS', 'DFP', 'SR1'. Default: 'BFGS'.
%
% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2006 by Miguel A. Carreira-Perpinan

function [X,F] = Obfgs(f,paramf,x0,tol,maxit,qN,H0)

% ---------- Argument defaults -----------
if ~exist('tol','var') | isempty(tol) tol = 1e-5; end;
if ~exist('maxit','var') | isempty(maxit) maxit = 10000; end;
n = length(x0);
if ~exist('H0','var') | isempty(H0) H0 = eye(n,n); end;
if ~exist('qN','var') | isempty(qN) qN = 'BFGS'; end;
% ---------- End of "argument defaults" ----------


x = x0;	H = H0;			% Current iterate
[ff g] = f(x',paramf{:});	% Current value of function and gradient
X = x'; F = ff;			% Lists of iterates and their function values
k = 0;				% Number of iterations
ff_old = 1e+5;                  % add by cqin
x_old = zeros(size(x));         % add by cqin

%fprintf(1,'Cycle %d: ff=%f  norm(g)=%f \n',k,ff,norm(g));

cnt = 0; max_cnt = 10;

%while (norm(g) > tol) & (k < maxit)  % original
%while (norm(g) > tol) & (norm(x-x_old) > 1e-7) & (k < maxit)
%while (abs(ff-ff_old) > 1e-6*abs(ff_old)) & (norm(x-x_old) > 1e-5) & (k < maxit)
%while (abs(ff-ff_old) > 1e-6) & (norm(x-x_old) > 1e-6) & (k < maxit)
while (abs(ff-ff_old) > tol*abs(ff_old)) & (k < maxit)
  
  ff_old = ff;                  % add by cqin
  x_old = x;                    % add by cqin
    
  p = -H*g';			% New search direction
  
  % Backtracking line search
  alpha = linesearch(f,paramf,x,p,ff,g);

  s = alpha*p;
  x = x + s;			% New iterate
  oldg = g;
  [ff g] = f(x',paramf{:});	% New value of function and gradient
  y = g - oldg;
  
  switch qN
   case 'DFP'
    Hy = H*y'; H = H - (Hy*Hy')/(y*Hy) + (s*s')/(y*s);
   case 'BFGS'
    r = y*s;
    if r > 0 H = (eye(n,n)-(s*y)/r)*H*(eye(n,n)-(y'*s')/r) + (s*s')/r; end;
   case 'SR1'
    % SR1 often generates indefinite H, so it should be implemented with
    % a trust-region strategy
    sHy = s - H*y';
    if abs(y*sHy) >= norm(y)*norm(sHy)*1e-8 H = H + (sHy*sHy')/(y*sHy); end;
  end
    
  %if norm(g)==norm(oldg)  % original
  if (abs(norm(g)-norm(oldg))<1e-6)
    cnt = cnt + 1;
    if (cnt>max_cnt) break; end;
  else
    cnt = 0;    % reset
  end
    
  X = [X;x']; F = [F;ff];	% Update lists
  k = k + 1;
    
  %fprintf(1,'Cycle %d: ff=%f norm(g)=%f norm(dx)=%f alpha=%f r=%f\n',k,ff,norm(g),norm(x-x_old),alpha,r);
  
end

if k==maxit
  fprintf(1,'Max iterations reached!\n');
end
