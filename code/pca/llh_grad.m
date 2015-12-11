% [a,s,grad] = llh_grad(S,N,lambda,phi) Factor analysis log-likelihood gradient
%
% Computes the log-likelihood gradient and two norms of it under a multinormal
% (Gaussian) model of mean 0 and covariance matrix lambda*lambda'+diag(phi)
% for a data set of N points with sample covariance S.
%
% In:
%   S: DxD covariance matrix of the data. S = TT'/N - E{T}E{T}' if T is
%      a DxN matrix with N data vectors as column vectors.
%   N: number of data points.
%   lambda: DxL matrix with factor loadings for L-factor solution.
%   phi: D-vector with uniquenesses for L-factor solution.
% Out:
%   a: maximum absolute value of the gradient.
%   s: square root of the average squared value of the gradient.
%   grad: Dx(L+1) matrix with the gradient of the log-likelihood.
%
% The gradient is (see Morrison: Multivariate Statistical Methods, 1990)
% the Dx(L+1) matrix:
%    -N [ Sigma_i*(I-S*Sigma_i)*lambda 0.5*diag(Sigma_i(I-S*Sigma_i)) ]
% where Sigma_i is the inverse of the reconstructed covariance matrix
% Sigma = lambda*lambda'+diag(phi).
%
% See also llh_normal.

% Copyright (c) 1997 by Miguel A. Carreira-Perpinan

function [a,s,grad] = llh_grad(S,N,lambda,phi)

[D,L] = size(lambda);

% Inverse reconstructed covariance matrix Sigma_i
if any(phi<=0)	% Then let's try to invert DxD matrix
  Sigma_i = inv( lambda*lambda'+diag(phi) );
else			% Otherwise just invert LxL matrix
  phi_i = diag(1./phi);
  Delta = inv( eye(L) + lambda' * phi_i * lambda );
  delta = phi_i * lambda * Delta;
  Sigma_i = ( eye(D) - delta * lambda' ) * phi_i;
end

% Log-likelihood gradient
grad = -N * [ Sigma_i*(eye(D)-S*Sigma_i)*lambda ...
             0.5*diag(Sigma_i*(eye(D)-S*Sigma_i)) ];

% Compute norms of the log-likelihood gradient
a = max(max(abs(grad)));
s = sqrt(trace(grad'*grad)/(D*(L+1)));
