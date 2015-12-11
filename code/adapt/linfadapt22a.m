% [Cy,dy,E] = linfadapt22a(X,Y,f,Ax,bx) Adaptation of a linear predictive mapping
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

function [Cy,dy,E] = linfadapt22a(X,Y,f,Ax,bx)

N = size(X,1); K = size(X,2)/2; P = size(Y,2)/2;

% Unpack arguments
switch f.type
 case 'rbf'
  MU = f.C; s = f.s; W = f.W(:,1:end-1)'; w = f.W(:,end)';
  % Transpose matrices, to conform to math notation
  MU = MU'; W = W'; w = w'; 
 case 'linf'
  W = f.W; w = f.w;
end

% Transpose matrices, to conform to math notation
X = X'; Y = Y'; 

PPy = zeros(2*P,4*P);
PP2y = zeros(4*P,4*P); 

tmp1 = zeros(4*P,1);
tmp2 = zeros(2*P,1);
Q = -eye(2*P);

for n=1:N
  
  % Linearly transformed Y and X
  Xn = X(:,n); Xn2 = reshape(Xn,2,K);
  Yn = Y(:,n); Yn2 = reshape(Yn,2,P);
  
  blkx = zeros(2*K,4*K); 
  for k=1:K
    %blkx(2*k-1:2*k,4*(k-1)+1:4*k) = kron(Xn2(:,k)',eye(2)); % original version
    blkx(2*k-1:2*k,4*(k-1)+1:4*k) = [Xn2(1,k) 0 Xn2(2,k) 0;0 Xn2(1,k) 0 Xn2(2,k)]; % fast version
  end
  
  XX = blkx * Ax' + bx'; % g(X)
  switch f.type
    case 'rbf'
     Phi = exp(-sqdist(MU',XX')/(2*s*s));	% Mx1, phi(g(Xn)) is column n
     Ytn = W * Phi + w; % f(g(X))  
   case 'linf'
     Ytn = W * XX + w; % f(g(X))
  end  
  Ytn2 = reshape(Ytn,2,P);
  
  blky = zeros(2*P,4*P); 
  for p=1:P
    %blky(2*p-1:2*p,4*(p-1)+1:4*p) = kron(Yn2(:,p)',eye(2));
    blky(2*p-1:2*p,4*(p-1)+1:4*p) = [Ytn2(1,p) 0 Ytn2(2,p) 0;0 Ytn2(1,p) 0 Ytn2(2,p)]; % fast version
  end
  
  Pny = -blky;
  PPy = PPy + Pny;
  PP2y = PP2y + Pny'*Pny;
  tmp1 = tmp1 + Pny'*Yn;
  tmp2 = tmp2 + Yn;
end

tmp = [PP2y PPy'*Q; Q'*PPy N*(Q'*Q)];
Cd = tmp \ ([-tmp1; -Q'*tmp2]);
Cy = (Cd(1:4*P))'; dy = (Cd(4*P+1:6*P))'; 

Ab = [Ax bx Cy dy]';

switch f.type
 case 'rbf'
  E = Ferr2_rbf2(Ab',f,X',Y',0);
 case 'linf'
  E = Ferr2_linf2(Ab',f,X',Y',0);
end

