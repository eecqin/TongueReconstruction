% [A,b,E] = linfadapt(X,Y,f,l) Adaptation of a linear predictive mapping
%
% Optimal linear transformation (given by A (2x2) and b (2x1)) wrt the
% reconstruction error E, for a linear predictive mapping f:
% - Linear mapping y = f(x) = W.x+w
% - 2D point transformation: xx = g(x) = A.x+b, yy = g(y) = A.y+b
% - Reconstruction error E(A,b) = \sum_n{|g(yn) - f(g(xn))|²}
%
% In:
%   X: Nx(2.K) matrix containing N 2K-dim contour points rowwise.
%   Y: Nx(2.P) matrix containing N 2P-dim contour points rowwise.
%   f: struct with fields W, w for the linear mapping (see linftrain.m).
%   l: regularisation parameter to minimise cond(A). Default: 0.
% Out:
%   A, b: optimal A (2x2 matrix) and b (2x1 vector).
%   E: error after adaptation.

% Copyright (c) 2009 by Miguel A. Carreira-Perpinan

function [A,b,E] = linfadapt(X,Y,f,l)

% Unpack arguments
W = f.W; w = f.w;

N = size(X,1); K = size(X,2)/2; P = size(Y,2)/2;

% Transpose matrices, to conform to math notation
X = X'; Y = Y';

PP = zeros(2*P,4); PP2 = zeros(4,4);
Q = kron(ones(P,1),eye(2)); for i=1:K Q = Q - W(:,2*i-1:2*i); end;

for n=1:N
  
  Xn = X(:,n); Xn2 = reshape(Xn,2,K);
  Yn = Y(:,n); Yn2 = reshape(Yn,2,P);
  
  PPn = kron(Yn2',eye(2));
  for i=1:K PPn = PPn - kron(Xn2(:,i)',W(:,2*i-1:2*i)); end;

  PP = PP + PPn;
  PP2 = PP2 + PPn'*PPn;
  
end

if exist('l','var') & ~isempty(l)
  % Test two cases for sign(det(A))
  PP2p = PP2 + l*(eye(4)-[0 0 0 1;0 0 -1 0;0 -1 0 0;1 0 0 0]);
  PP2n = PP2 + l*(eye(4)+[0 0 0 1;0 0 -1 0;0 -1 0 0;1 0 0 0]);
  tmpp = [PP2p PP'*Q;Q'*PP N*Q'*Q] \ ([PP N*Q]'*w);
  tmpn = [PP2n PP'*Q;Q'*PP N*Q'*Q] \ ([PP N*Q]'*w);
  if Ferr_linf(tmpp',f,X',Y',l) > Ferr_linf(tmpn',f,X',Y',l)
    tmp = tmpn;
  else
    tmp = tmpp;
  end
else
  tmp = [PP2 PP'*Q;Q'*PP N*Q'*Q] \ ([PP N*Q]'*w);
end
A = reshape(tmp(1:4),2,2); b = tmp(5:6);

E = Ferr_linf(tmp',f,X',Y',l);

