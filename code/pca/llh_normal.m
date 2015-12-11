% l = llh_normal(Sdat,N,Smod) Log-likelihood for normal model N(0,cov)
%
% Computes the log-likelihood under a multinormal (Gaussian) model of
% mean 0 and covariance matrix Smod for a data set of N points with sample
% covariance Sdat. If Smod is omitted, it gives llh_normal(Sdat,N,Sdat),
% which is equal to max{llh_normal(Sdat,N,S)} over all S.
%
% Sdat: DxD covariance matrix of the data. Sdat = TT'/N - E{T}E{T}' if T is
%       a DxN matrix with N data vectors as column vectors.
% N: number of data points.
% Smod: DxD covariance matrix of the model (default Sdat).
%
% By definition llh(Smod) = \sum^N_{n=1}{\ln{p(t_n|Smod)}}.
% For p ~ N(0,Smod) this simplifies to
% llh(Smod) = -N/2 { D\ln{2\pi} + \ln{|Smod|} + trace{Smod^{-1} Sdat}
% }.
%
% Observe that principal component analysis (PCA) and factor analysis
% (FA) fall into this category of models.
%
% See also llh_grad.

% Copyright (c) 1997 by Miguel A. Carreira-Perpinan

function l = llh_normal(Sdat,N,Smod)

% Argument defaults
if nargin==2 Smod=Sdat; end;

l = -N/2 * ( length(Sdat(1,:))*log(2*pi) + log(det(Smod)) + ...
    trace(inv(Smod)*Sdat) );
