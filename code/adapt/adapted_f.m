% [Y,fgX,gX] = adapted_f(X,f,A,b) Apply function f adapted with 2D transformation A,b
%
% Computes Y = invg(f(g(X))) where g(x) = A.x+b is applied to 2D blocks in X.
% A,b will have been obtained with:
% - linfadapt if using a linear predictive mapping f
% - rbfadapt if using an RBF predictive mapping f.
%
% In:
%   X: Nx(2.K) matrix containing N 2K-dim contour points rowwise.
%   f: struct for the function (see linftrain.m, rbftrain.m).
%   A, b: A (2x2 matrix) and b (2x1 vector) for the 2D transformation.
% Out:
%   Y: Nx(2.P) matrix containing N 2P-dim contour points rowwise.
%   fgX: Nx(2.P) matrix containing N 2P-dim predicted contour points.
%   gX: Nx(2.K) matrix containing N 2K-dim transformed contour points.

% Copyright (c) 2009 by Miguel A. Carreira-Perpinan

function [Y,fgX,gX] = adapted_f(X,f,A,b)

K = size(X,2)/2;
fcn = str2func(f.type);			% Use string as function handle

gX = bsxfun(@plus,kron(eye(K),A)*X',kron(ones(K,1),b));		% g(X)
fgX = fcn(gX',f);						% f(g(X))
P = size(fgX,2)/2;
Y = bsxfun(@minus,fgX,kron(ones(P,1),b)') / kron(eye(P),A)';	% invg(f(g(X)))




% The matrix multiplication and linear system involving kron() could be
% implemented more efficiently with a loop over the 2D blocks, but for
% small K or P this is good enough.

