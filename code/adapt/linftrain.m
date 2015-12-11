% [f,fX,E] = linftrain(X,Y) Train linear function y = f(x) = W.x+w
%
% In:
%   X: NxL matrix, N L-dim data points rowwise.
%   Y: NxD matrix, N D-dim data points rowwise.
% Out:
%   f: (struct) the linear function, with fields:
%      type='linf', W (DxL), w (Dx1).
%   fX: NxD matrix, f(X).
%   E: fit error.

% Copyright (c) 2009 by Miguel A. Carreira-Perpinan

function [f,fX,E] = linftrain(X,Y)

N = size(X,1);
f.type = 'linf';

X1 = sum(X,1)'; XX = X'*X - X1*(X1'/N);
f.W = (Y'*X-sum(Y,1)'*X1'/N) / XX; f.w = sum(Y-X*f.W',1)'/N;

if nargout>1 fX = linf(X,f); end
if nargout>2 E = sum(sum((Y-fX).^2)); end

