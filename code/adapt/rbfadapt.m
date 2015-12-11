% [A,b,E] = rbfadapt(X,Y,f,A,b,l,tol) Adaptation of an RBF predictive mapping
%
% Optimal linear transformation (given by A (2x2) and b (2x1)) wrt the
% reconstruction error E, for an RBF predictive mapping f:
% - RBF mapping y = f(x) = W.phi(x)+v, phi(x)_m = exp(-|(x-mu_m)/s|²/2)
% - 2D point transformation: xx = g(x) = A.x+b, yy = g(y) = A.y+b
% - Reconstruction error E(A,b) = \sum_n{|g(yn) - f(g(xn))|²}
%
% In:
%   X: Nx(2.K) matrix containing N 2K-dim contour points rowwise.
%   Y: Nx(2.P) matrix containing N 2P-dim contour points rowwise.
%   f: struct with fields C, s, W for the RBF (see rbftrain.m).
%   A, b: initial A (2x2 matrix) and b (2x1 vector).
%   l: regularisation parameter to minimise cond(A). Default: 0.
%   tol: convergence tolerance for BFGS. Default: 1e-5.
% Out:
%   A, b: optimal A (2x2 matrix) and b (2x1 vector).
%   E: error after adaptation.

% Copyright (c) 2009 by Miguel A. Carreira-Perpinan

function [A,b,E] = rbfadapt(X,Y,f,A,b,l,tol)

% ---------- Argument defaults ----------
if ~exist('A','var') | isempty(A) | ~exist('b','var') | isempty(b)
  A = eye(2); b = [0;0];	% Initial A, b
end;
if ~exist('l','var') l = []; end;
if ~exist('tol','var') | isempty(tol) tol = 1e-5; end;
% ---------- End of "argument defaults" ----------

[tmp1,tmp2,H0] = Ferr_rbf([A(:);b(:)]',f,X,Y,l); % Initial Hessian estimate
[tmp1,tmp2] = Obfgs(@Ferr_rbf,{f,X,Y,l},[A(:);b(:)],tol,[],[]);
A = reshape(tmp1(end,1:4),2,2); b = tmp1(end,5:6)';
E = tmp2(end);

