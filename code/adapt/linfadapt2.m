% [A,b,E] = linfadapt2(X,Y,f,l) Adaptation of a linear predictive mapping
%
% Optimal linear transformation (given by A (2x2) and b (2x1)) wrt the
% reconstruction error E, for a linear predictive mapping f:
% - Linear mapping y = f(x) = W.x+w
% - 2D point transformation: xx = g(x) = A.x+b, yy = g(y) = A.y+b
% - Reconstruction error E(A,b) = \sum_n{|g(yn) - f(g(xn))|�}
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

function [Ab,E] = linfadapt2(X,Y,f,l)

% Unpack arguments
W = f.W; w = f.w;

N = size(X,1); K = size(X,2)/2; P = size(Y,2)/2;

% Transpose matrices, to conform to math notation
X = X'; Y = Y';

PPx = zeros(2*P,4*K); PPy = zeros(2*P,4*P);
PP2x = zeros(4*K,4*K); PP2y = zeros(4*P,4*P); PP2xy = zeros(4*K,4*P);

for n=1:N
  
  % Linearly transformed Y and X
  Xn = X(:,n); Xn2 = reshape(Xn,2,K);
  Yn = Y(:,n); Yn2 = reshape(Yn,2,P);
  blkx = zeros(2*K,4*K); blky = zeros(2*P,4*P); 
  
  for k=1:K
    blkx(2*k-1:2*k,4*(k-1)+1:4*k) = kron(Xn2(:,k)',eye(2));
  end

  for p=1:P
    blky(2*p-1:2*p,4*(p-1)+1:4*p) = kron(Yn2(:,p)',eye(2));
  end
  
  Pnx = -W * blkx; Qx = -W;
  Pny = blky; Qy = eye(2*P);
  
  PPx = PPx + Pnx; PPy = PPy + Pny;  
  PP2x = PP2x + Pnx'*Pnx; 
  PP2y = PP2y + Pny'*Pny;
  PP2xy = PP2xy + Pnx'*Pny;
  
end

tmp = [PP2x PPx'*Qx PP2xy PPx'*Qy; ...
       Qx'*PPx N*(Qx'*Qx) Qx'*PPy N*Qx'*Qy; ...
       PP2xy' PPy'*Qx PP2y PPy'*Qy; ...
       Qy'*PPx N*Qy'*Qx Qy'*PPy N*Qy'*Qy];
Ab = tmp \ ([PPx N*Qx PPy N*Qy]'*w);

% $$$ Ax = (Ab(1:4*K))'; bx = (Ab(4*K+1:6*K))'; 
% $$$ Ay = (Ab(6*K+1:6*K+4*P))'; by = (Ab(6*K+4*P+1:end))';

Ab = Ab';
E = Ferr_linf2(Ab,f,X',Y');

