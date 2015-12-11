% [A,b] = Flin(X,Y,W,v) Adaptation (linear prediction): optimal transformation
%
% Optimal linear transformation (given by A (2x2) and b (2x1)) wrt the
% reconstruction error F, for a linear prediction mapping:
% - Linear mapping y = f(x) = W.x+v
% - 2D point transformation: xx = g(x) = A.x+b, yy = g(y) = A.y+b
% - Reconstruction error F(A,b) = \sum_n{|g(yn) - f(g(xn))|²}
%
% In:
%   X: Nx(2.K) matrix containing N 2K-dim contour points rowwise.
%   Y: Nx(2.P) matrix containing N 2P-dim contour points rowwise.
%   W: (2.K)x(2.P) matrix containing 2K 2P-dim weight vectors rowwise.
%   v: 1x(2.P) vector containing the 2P-dim bias vector rowwise.
% Out:
%   A, b: optimal A (2x2 matrix) and b (2x1 vector).

% Copyright (c) 2008 by Miguel A. Carreira-Perpinan

function [A,b] = Flin(X,Y,W,v)

N = size(X,1); K = size(X,2)/2; P = size(Y,2)/2;

% Transpose matrices, to conform to math notation
X = X'; Y = Y'; W = W'; v = v';

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

tmp = [PP2 PP'*Q;Q'*PP N*Q'*Q] \ [PP N*Q]'*v;
A = reshape(tmp(1:4),2,2); b = tmp(5:6);

