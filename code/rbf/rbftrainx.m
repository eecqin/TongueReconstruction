% [f,fXv,E,Er] = rbftrainx(X,Y,Xv,Yv,C,s,l) Crossvalidate RBF to map f(X)=Y
%
% Train on (X,Y) crossvalidating the width and regularisation parameter on
% (Xv,Yv).
%
% In:
%   X: NxL matrix, N L-dim data points rowwise.
%   Y: NxD matrix, N D-dim data points rowwise.
%   Xv: NNxL matrix, NN L-dim data points rowwise.
%   Yv: NNxD matrix, NN D-dim data points rowwise.
%   C: if a scalar, the #BFs to use;
%      else, MxL matrix, M L-dim initial values of the RBF centres rowwise.
%   s: list of SS RBF widths.
%   l: list of LL regularisation parameter values.
% Out:
%   f: (struct) the RBF with smallest fit error E on (Xv,Yv), with fields:
%      type='rbf', centres C (MxL), width s, weights W (Dx(M+1), incl. bias),
%      regularisation parameter l.
%   fXv: NxD matrix, f(Xv).
%   E,Er: SSxLL arrays, fit & regularisation errors on (Xv,Yv).

% Copyright (c) 2009 by Miguel A. Carreira-Perpinan

function [f,fXv,E,Er] = rbftrainx(X,Y,Xv,Yv,C,s,l)

N = size(X,1); ff.type = 'rbf';

% Fit centres
if isscalar(C)
  M = C; e = Inf;	% run kmeans 10 times with rnd init and pick best
  for i=1:10
    [tmpC,tmp,tmpe] = kmeans(X,M);
    if tmpe(end) < e C = tmpC; e = tmpe(end); end
  end
else
  C = kmeans(X,size(C,1),C);
end
ff.C = C; sqdCX = sqdist(C,X); M = size(C,1);

% Fit an RBF for each combination of (s,l) given the centres
Eold = Inf; E = zeros(length(s),length(l)); Er = E;
warning('off','MATLAB:nearlySingularMatrix');
for i = 1:length(s)
  % Amortise computation of Gram matrices for different l and fixed s
  G = exp(-sqdCX/(2*s(i)*s(i))); G1 = sum(G,2); GG = G*G' - G1*(G1'/N);
  YG = Y'*G'-mean(Y',2)*G1';
  for j = 1:length(l)
    W = YG / (GG+spdiags(repmat(l(j),M,1),0,M,M)); w = mean(Y'-W*G,2);
    ff.s = s(i); ff.W = [W w]; ff.l = l(j); ffXv = rbf(Xv,ff);
    E(i,j) = sum(sum((Yv-ffXv).^2)); Er(i,j) = l(j)*(W(:)'*W(:));
    if E(i,j) < Eold f = ff; fXv = ffXv; Eold = E(i,j); end
  end
end
warning('on','MATLAB:nearlySingularMatrix');

