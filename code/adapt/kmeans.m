% [mu,labels,errfunc,code] = kmeans(X,K[,init,maxit,tol])
% K-means clustering of dataset X
%
% In:
%   X: NxD matrix containing N D-dimensional data points rowwise.
%   K: integer in [1,N] containing the desired number of clusters.
%   init: initialisation, one of 'rndlabels' (random assignment),
%      'rndmeans' (random means in the range of X), 'rndsubset' (random
%      subset of data points) or a KxD matrix containing the initial K
%      cluster means. Default: 'rndmeans'.
%   maxit: maximal number of iterations. Default: Inf.
%   tol: small positive number, tolerance in the relative change of the
%      error function to stop iterating. Default: 0.
% Out:
%   mu: KxD matrix containing the K cluster means.
%   labels: Nx1 list containing the cluster labels (1 to K). Thus,
%      labels(n) = k says point Xn is in cluster k.
%   errfunc: list of values (for each iteration) of the error function.
%      This is the sum over all clusters of the within-cluster sums of
%      point-to-mean distances.
%   code: stopping code: 0 (tolerance achieved), 1 (maxit exceeded).
%
% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2003 by Miguel A. Carreira-Perpinan

function [mu,labels,errfunc,code] = kmeans(X,K,init,maxit,tol)

% ---------- Argument defaults ----------
if ~exist('init','var') | isempty(init) init = 'rndmeans'; end;
if ~exist('maxit','var') | isempty(maxit)
% $$$   maxit = max(100,floor(size(X,1)/100));
  maxit = Inf;
end;
if ~exist('tol','var') | isempty(tol) tol = 0; end;
% ---------- End of "argument defaults" ----------

% Initialisation
mu = km_init(X,K,init);
[labels,errfunc] = km_le(X,mu);

% Iteration
code = -1; it = 0;
while code < 0
  oldlabels = labels; oldmu = mu;		% Keep if back up needed
  mu = km_m(X,labels,mu);			% New means,
  [labels,errfunc(end+1)] = km_le(X,mu);	% labels and error value
  it = it + 1;
  % Stopping condition
  if errfunc(end) >= (1-tol)*errfunc(end-1)
    code = 0;
  elseif it >= maxit
    code = 1;
  end
end

% Back up if the last iteration increased the error function. This can
% happen due to numerical inaccuracy if two points Xn are very close.
if errfunc(end) > errfunc(end-1)
  it = it - 1;
  errfunc = errfunc(1:end-1);
  mu = oldmu;
  labels = oldlabels;
end


% [labels,e] = km_le(X,mu) Labels and error function given means.
%
% In:
%   X, mu: see above.
% Out:
%   labels: see above.
%   e: value of the error function.

function [labels,e] = km_le(X,mu)

sqd = sqdist(X,mu);
[e,labels] = min(sqd,[],2);	% Assignment of data points to means
e = sum(e);			% Error


% mu = km_m(X,labels,mu[,deadmu]) Means given labels.
%
% In:
%   X, labels, mu: see above.
%   deadmu: what to do with a dead mean (which has no associated data
%      points), one of 'rnd' (assign to it a random data point) or
%      'leave' (leave as is). Default: 'rnd'.
% Out:
%   mu: new means.

function mu = km_m(X,labels,mu,deadmu)

% Argument defaults
if ~exist('deadmu','var') | isempty(deadmu) deadmu = 'rnd'; end;

for k=1:size(mu,1)
  tmp = find(labels==k); ltmp = length(tmp);
  if ltmp > 0
    mu(k,:) = mean(X(tmp,:),1);
  elseif strcmp(deadmu,'rnd')
    % Dead mean: assign to it one data point at random
    mu(k,:) = X(floor(rand*size(X,1)+1),:);
    % Otherwise:
    % Leave dead mean as is (works badly because often means end up dead)
  end
end


% mu = km_init(X,K,init)
%
% Initialise means.
%
% In:
%   X, K, init: as above.
% Out:
%   mu: as above.

function mu = km_init(X,K,init)

[N,D] = size(X);
if ischar(init)
%  rand('state',sum(100*clock));
  switch init
   case 'rndlabels'
    % labels = random assignment
    mu = km_m(X,[(1:K)';floor(K*rand(N-K,1)+1)],zeros(K,D));
   case 'rndmeans'
    % mu = uniform random points in the range of X
    m = min(X,[],1); M = max(X,[],1);
    mu = bsxfun(@plus,m,bsxfun(@times,M-m,rand(K,D)));
   case 'rndsubset'
    % mu = random subset of data points (works badly because often means
    % compete inside clusters)
    tmp = randperm(N); mu = X(tmp(1:K),:);
  end
else
  % mu = user-provided
  mu = init;
end

