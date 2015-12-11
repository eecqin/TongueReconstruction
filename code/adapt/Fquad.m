% function [f,g,H] = Fquad(X,A,b,c) Quadratic function
%
% f(x) = ½x'*A*x + b'*x + c
%
% In:
%   X: N x n list of row vectors.
%   A: n x n symmetric matrix.
%   b: n x 1 vector.
%   c: real number.
% Out:
%   f: N x 1 list of function values.
%   g: N x n list of gradient vectors.
%   H: n x n x N list of Hessian matrices.

% Copyright (c) 2005 by Miguel A. Carreira-Perpinan

function [f,g,H] = Fquad(X,A,b,c)

N = size(X,1);

f = sum(X.*(X*A')/2,2) + X*b + c;

if nargout > 1
  g = X*A'+repmat(b',N,1);
  if nargout > 2
    H = repmat(A,[1 1 N]);
  end
end

