% alpha = linesearch(f,paramf,x,p[,ff,g,alpha0,rho,c]) Backtracking line search
%
% Reference: procedure 3.1, p. 41ff in:
%   Nocedal and Wright: "Numerical Optimization", Springer, 1999.
%
% This procedure is intended to be called from an optimisation procedure
% such as Osteepdesc.
%
% In:
%   f, paramf: see Osteepdesc.m.
%   x: n x 1 vector, the current iterate.
%   p: n x 1 vector, the search direction.
%   ff, g: value at x of the function and gradient as returned by f.
%      Default: compute them.
%   alpha0: initial step size. Default: 1.
%   rho: rate of decrease of the step size. Default: 0.8.
%   c: Wolfe condition. Default: 1e-4.
% Out:
%   alpha: step size.
%
% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2006 by Miguel A. Carreira-Perpinan

function alpha = linesearch(f,paramf,x,p,ff,g,alpha0,rho,c)

% ---------- Argument defaults ----------
if ~exist('ff','var') | isempty(ff) | ~exist('g','var') | isempty(g)
  [ff g] = f(x',paramf{:});
end;
if ~exist('alpha0','var') | isempty(alpha0) alpha0 = 1; end;
if ~exist('rho','var') | isempty(rho) rho = 0.8; end;
if ~exist('c','var') | isempty(c) c = 1e-4; end;
% ---------- End of "argument defaults" ----------

alpha = alpha0;
tmp = c*g*p;
while f((x+alpha*p)',paramf{:}) > ff + alpha*tmp
  alpha = rho*alpha;
end

