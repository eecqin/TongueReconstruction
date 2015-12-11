% [U,V,llh,s2,Sigma] = pca2(S,N,L) Principal Component Analysis of cov. matrix
%
% Principal component analysis of a symmetric matrix S is the singular 
% value decomposition S = U.V.U' with U orthogonal and V diagonal. To
% extract only L components, U and V are restricted accordingly.
%
% In:
%   S: DxD covariance matrix (>0).
%   N: sample size.
%   L: number of principal components.
% Out:
%   U: DxD matrix containing the eigenvectors of S by columns.
%   V: Dx1 vector containing the associated eigenvalues of S.
%   llh: log-likelihood of the probabilistic PCA model.
%   s2: uniqueness parameter that optimises the log-likelihood of the
%      probabilistic PCA model.
%   Sigma: DxD symmetric matrix containing the normal model covariance
%      matrix for PCA, i.e. lambda*lambda'+I/s2. This is useful to
%      compute the log-likelihood under the same model of a different
%      data set using the function llh_normal.
%
% The probabilistic PCA model is described in:
%
%   M. E. Tipping and C. M. Bishop: "Probabilistic Principal Component
%   Analysis", Tech. rep. NCRG/97/010, Aston University, 1997.
%
% See also pca, pca_test, fa, lin_proj_error, llh_normal.

% Copyright (c) 1997 by Miguel A. Carreira-Perpinan

function [U,V,llh,s2,Sigma] = pca2(S,N,L);

[D,D] = size(S);

% Eigenvectors and eigenvalues
[U,Vtemp,U] = svd(S);
U = U(:,1:L);
Vtemp = diag(Vtemp);
V = Vtemp(1:L);

if nargout > 2
  if L<D
    s2 = sum(Vtemp(L+1:D))/(D-L);
    % The log-likelihood could also be computed with
    % llh_normal(S,N,Sigma) but more slowly.
    llh = -N/2*(D*log(2*pi)+(D-L)*log(s2)+sum(log(V))+L+sum(Vtemp(L+1:D))/s2);
    %llh=-N/2*(D*(1+log(2*pi))+(D-L)*log(sum(Vtemp(L+1:D))/(D-L))+sum(log(V)));
    Sigma = s2*eye(D) + U*diag(V-s2)*U';
  else
    s2=0;
    llh=-N/2*(D*(1+log(2*pi))+sum(log(Vtemp)));
    Sigma=S;
  end
end
