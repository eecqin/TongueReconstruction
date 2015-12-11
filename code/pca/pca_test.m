% [reject,stat,chi,nu] = pca_test(eigenval,N,L) Test for PCA
%
% Computes a likelihood ratio test statistic and, using the
% chi-squared, tells whether the null hypothesis:
%	H0: the data can be explained by *exactly* L principal components
% should be rejected or not. The significance of the goodness-of-fit test
% is fixed to 95%.
%
% In:
%   eigenval: Dx1 vector containing the eigenvalues of the PCA.
%   N: number of data vectors (necessary to compute the log-likelihood).
%   L: number of principal components (>=0).
% Out:
%   reject: 1 if H0 is rejected, 0 otherwise. reject = (stat > chi). A -1
%       value is returned if the resulting number of degrees of freedom is
%       not positive.
%   stat: the likelihood ratio test statistic numerical value.
%   chi: the chi-squared critical point.
%   nu: number of degrees of freedom of the chi-squared.
%
% See also pca, pca2, chi2inv.

% Copyright (c) 1997 by Miguel A. Carreira-Perpinan

function [reject,stat,chi,nu] = pca_test(eigenval,N,L)

D = length(eigenval);

nu = (D-L+2)*(D-L-1)/2;	% Degrees of freedom of the chi-square

if nu<=0
  reject = -1;
else
  chi = chi2inv(0.95,nu);
  stat = (N-(2*D+11)/6) * (D-L) * log( sum(eigenval(L+1:D))/(D-L) / ...
      prod(eigenval(L+1:D))^(1/(D-L)) );
  reject = stat > chi;
end
