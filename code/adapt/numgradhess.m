% [g,H,z,ge,He] = numgradhess(f,paramf,X[,h,tol]) Numerical gradient & Hessian
%
% Evaluates numerically (with a finite-difference approximation based on
% values of the function f) the gradient and Hessian of f at the points X,
% and optionally compares them with the ones returned by f itself (useful to
% confirm the latter are correct).
%
% Example of use:
%   n=2; X=randn(3,n); A=randn(n,n); A=A+A'; b=randn(n,1); c=randn;
%   h = 0; [g,H,z] = numgradhess(@Fquad,{A,b,c},X,h)
% The function Fquad.m (available separately) shows the format in which the
% function value, gradient and Hessian must be returned. Note that it must
% return the values of the function, gradient and Hessian for a list of N
% input vectors, eg [f,g,H] = Fquad(X,A,b,c) where X is N x n, even if you
% call numgradhess with a single input vector.
%
% Important notes:
% - The perturbation size e can be neither too large (large truncation error)
%   nor too small (large roundoff error). For a forward-difference
%   approximation to the gradient, e ~ sqrt(eps/2) is advisable.
% - For a second forward-difference approximation to the Hessian (option h=0),
%   the roundoff error can be terrible, so it needs a larger e.
% - I use e = 100*sqrt(eps/2) for both gradient and Hessian as a tradeoff, but
%   the Hessian estimate can often be pretty bad. A much better estimate can
%   be obtained by applying a forward-difference to the true gradient
%   expression (option h=1).
% - I also multiply e times the largest-magnitude component in each X(i,:)
%   vector. If the components vary widely in magnitude this can give bad
%   results; using a different e for each component would be better but is
%   messier.
% 
% Reference: p. 195ff in:
%   Nocedal and Wright: "Numerical Optimization", 2nd ed, Springer, 2006.
%
% In:
%   f: a handle to a function that takes as input a list of vectors and
%      returns 3 lists: the function values, the gradients and the Hessians.
%   paramf: cell array containing the parameters for f. Use {} to pass no
%      parameters, or to use f's default parameters.
%   X: N x n list of row vectors.
%   h: 0 to estimate the Hessian using only function evaluations; 1 to
%      estimate it using evaluations of the true gradient provided by f.
%      Default: 0.
%   tol: tolerance in the relative error between the numerical and the given
%      derivatives. Default: 1e3*e with e defined below.
% Out:
%   g: N x n list of numerical gradient vectors.
%   H: n x n x N list of numerical Hessian matrices.
% Optional output arguments:
%   z: Nx2 binary array with entry=0 if the numerical gradient and Hessian
%      (columns 1 and 2, resp.) agree with the ones provided by f to the
%      tolerance specified, entry=1 otherwise.
%   ge, He (): like g, H but containing the errors g - giveng, H - givenH.
%
% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2010 by Miguel A. Carreira-Perpinan

function [g,H,z,ge,He] = numgradhess(f,paramf,X,h,tol)

e = 100*sqrt(eps/2);	% a reasonable perturbation size
% ---------- Argument defaults ----------
if ~exist('h','var') | isempty(h) h = 0; end;
if ~exist('tol','var') | isempty(tol) tol = 1e3*e; end;
% ---------- End of "argument defaults" ----------

[N,n] = size(X);
g = zeros(N,n);
E = speye(n);		% perturbation directions for the gradient
if nargout > 1		% perturbation directions for the Hessian
  Ei = repmat(E,n,1); Ej = reshape(Ei,n,n^2)'; Eij = Ei+Ej; clear Ei Ej;
  H = zeros(n,n,N);
  if nargout > 2
    z = zeros(N,2);
    if nargout > 3
      ge = g;
      if nargout > 4
        He = H;
      end
    end
  end
end

for i=1:N
  
  x = X(i,:);
  ee = e*noz(max(abs(x)));		% perturbations propto xm

  if nargout > 2
    [fx,gg,HH] = f(x,paramf{:});	% given values of gradient/Hessian
  elseif h~=0
    [fx,gg] = f(x,paramf{:});		% given values of gradient
  else
    fx = f(x,paramf{:});
  end

  F = f(bsxfun(@plus,ee*E,x),paramf{:});
  
  g(i,:) = (F - fx)'/ee;			% numerical gradient (fwd diff)
  
  if nargout > 2				% gradient relative error
    tmpg = g(i,:) - gg;
    z(i,1) = any(abs(tmpg) > tol*abs(noz(gg)));
    if nargout > 3
      ge(i,:) = tmpg;
    end
  end
  
  if nargout > 1				% numerical Hessian (fwd diff)
    if h==0					% ...using only f
      H(:,:,i) = (bsxfun(...
        @minus,bsxfun(@minus,reshape(...
          f(bsxfun(@plus,ee*Eij,x),paramf{:}),n,n),F),F')+fx)/(ee^2);
    else
      [fxe,gge] = f(bsxfun(@plus,ee*E,x),paramf{:});	% ...using the
      H(:,:,i) = bsxfun(@minus,gge,gg)/ee;		% ...given gradient
      H(:,:,i) = (H(:,:,i) + H(:,:,i)')/2;		% symmetrise
    end
    if nargout > 2				% Hessian relative error
      tmpH = H(:,:,i) - HH;
      z(i,2) = any(abs(tmpH(:)) > tol*abs(noz(HH(:))));
      if nargout > 4
        He(:,:,i) = tmpH;
      end
    end
  end
  
end


function y = noz(x)

y = x; y(y==0) = 1;

