% [f,fXv,E,Er] = linftrainx(X,Y,Xv,Yv,l) Crossvalidate linear function to map f(X)=Y
%
% Train on (X,Y) crossvalidating the regularisation parameter on (Xv,Yv).
%
% In:
%   X: NxL matrix, N L-dim data points rowwise.
%   Y: NxD matrix, N D-dim data points rowwise.
%   Xv: NNxL matrix, NN L-dim data points rowwise.
%   Yv: NNxD matrix, NN D-dim data points rowwise.
%   l: list of LL regularisation parameter values.
% Out:
%   f: (struct) the linear function, with fields:
%      type='linf', W (DxL), w (Dx1), regularisation parameter l.
%   fXv: NxD matrix, f(Xv).
%   E,Er: LLx1 arrays, fit & regularisation errors on (Xv,Yv).

% Copyright (c) 2010 by Chao Qin

function [f,fXv,E,Er] = linftrainx(X,Y,Xv,Yv,l)

N = size(X,1); ff.type = 'linf';

Eold = Inf; E = zeros(length(l),1); Er = E;
for j = 1:length(l)  
  ff = linftrain(X,Y,l(j)); ffXv = linf(Xv,ff);
  E(j) = sum(sum((Yv-ffXv).^2)); Er(j) = l(j)*(ff.W(:)'*ff.W(:));
  if E(j) < Eold f = ff; fXv = ffXv; Eold = E(j); end
end
