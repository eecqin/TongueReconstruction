% [Y,J] = rbf(X,f) Value and Jacobian wrt x of RBF f(x)
%
% See rbftrain.
%
% In:
%   X: NxL matrix, N L-dim data points rowwise.
%   f: (struct) the RBF.
% Out:
%   Y: NxD matrix, N D-dim outputs Y = f(X).
%   J: DxL Jacobian matrix (assumes N=1 input only).

% Copyright (c) 2010 by Miguel A. Carreira-Perpinan

function [Y,J] = rbf(X,f)

phi = exp(-sqdist(f.C,X)/(2*f.s*f.s));
Y = (f.W*[phi;ones(1,size(X,1))])';
if nargout>1
  J = (f.W(:,1:end-1)*(diag(sparse(phi))*bsxfun(@minus,f.C,X)))/(f.s^2);
end

