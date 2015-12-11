% [U,V,m,S,llh,Sigma,s2,error2] = pca(X,L) Principal Component Analysis
%
% Computes the L principal components of the data matrix X.
%
% In:
%   X: NxD data matrix of N D-dimensional row vectors.
%   L: number of principal components.
% Out:
%   U: DxD matrix containing the principal components of X by columns.
%   V: Dx1 vector containing the associated eigenvalues.
%   m: Dx1 vector with the sample mean.
%   S: DxD sample covariance matrix.
%   llh: llh(i) is the log-likelihood of the probabilistic PCA model
%      using i components.
%   Sigma: DxD symmetric matrix containing the normal model covariance
%      matrix for PCA, i.e. lambda*lambda'+I/s2. This is useful to
%      compute the log-likelihood under the same model of a different
%      data set using the function llh_normal.
%   s2: s2(i) is the uniqueness parameter that optimises the
%      log-likelihood of the probabilistic PCA model using i
%      components.
%   error2: error2(i) is the quadratic average sample error using i
%   components.
%
% The probabilistic PCA model is described in:
%
%   M. E. Tipping and C. M. Bishop: "Probabilistic Principal Component
%   Analysis", Tech. rep. NCRG/97/010, Aston University, 1997.
%
% See also pca2, pca_test, fa, lin_proj_error, llh_normal.

% Copyright (c) 1997 by Miguel A. Carreira-Perpinan

function [U,V,m,S,llh,Sigma,s2,error2] = pca(X,L);

[N,D] = size(X);
m = mean(X);                      % Sample mean
X = X - ones(N,1)*m;              % Centre X
S = X'*X/N;                       % Sample covariance
[U,V,llh,s2,Sigma] = pca2(S,N,L);

if nargout>7
  error2 = lin_proj_error(X,U);
end
