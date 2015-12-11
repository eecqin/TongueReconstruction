% [e2,Xrec,eoo] = lin_proj_error(X,U,V) Sample average error with linear projection
%
% Given a set of N D-dimensional vectors X(NxD), a matrix U(DxL) and a 
% matrix V(DxL), the reconstructed set of vectors by projecting with U
% and backprojecting with V is XUV'. lin_proj_error computes the
% average quadratic error incurred, i.e. norm{X-XUV'}^2_2 / N, and the
% average "maximum norm" error, i.e. norm{X-XUV'}^2_oo / N.
%
% X is assumed to be zero-centred.
%
% In:
%   X: NxD data matrix of N D-dimensional row vectors.
%   U: DxL matrix (projection to L-space).
%   V: DxL matrix (backprojection to D-space). Defaults to the
%      transposed pseudoinverse of U.
% Out:
%   e2: the average quadratic error for the sample X.
%   Xrec: NxD data matrix containing the reconstructed data set.
%   eoo: the average "maximum norm" error for the sample X.
%
% See also pca, fa.

% Copyright (c) 1998 by Miguel A. Carreira-Perpinan

function [e2,Xrec,eoo] = lin_proj_error(X,U,V)

% Argument defaults
if nargin==2 V=pinv(U)'; end;

[N,D] = size(X);
Xrec = X*U*V';
% Reconstruction mean squared error
% This is much faster than doing e = sum(diag((X-Xrec)'*(X-Xrec)))/N;
e2 = sum(sum((X-Xrec).^2))/N;
if nargout > 2
  % Reconstruction mean "maximum norm" error
  eoo = sum(max((X-Xrec)'))/N;
end
