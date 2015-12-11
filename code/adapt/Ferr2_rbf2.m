% [ff,gg,HH] = Ferr2_rbf2(Ab,f,X,Y,l) Local adaptation error with RBF mapping
%
% Reconstruction error E(A,b,C,d) wrt the linear transformation (given by A (2.Kx2.K),
% b (2.Kx1), C (2.Px2.P), d (2.Px1)), for an RBF predictive mapping:
% - RBF mapping y = f(x) = W.phi(x)+w, phi(x)_m = exp(-|(x-mu_m)/s|²/2)
% - 2D point transformation: xx = gx(x) = A.x+b, yy = gy(y) = C.y+d
% - Reconstruction error E(A,b) = \sum_n{|g(yn) - f(g(xn))|²}
%
% In:
%   Ab: 1x(6.K+6.P) vector
%   f: struct with fields C, s, W for the RBF (see rbftrain.m).
%   X: Nx(2.K) matrix containing N 2K-dim contour points rowwise.
%   Y: Nx(2.P) matrix containing N 2P-dim contour points rowwise.
%   l: regularisation parameter to minimise cond(A). Default: 0.
% Out:
%   ff: error value.
%   gg: 1x(6.K+6.P) error gradient.
%   HH: (6.K+6.P)x(6.K+6.P) error approximate Hessian (a la Gauss-Newton).

% Copyright (c) 2010 by Chao Qin and Mohsen Farhadloo

function [ff,gg,HH] = Ferr2_rbf2(Ab,f,X2,Y2,l)
  
NN = size(Ab,1);
ff = []; gg = []; HH = [];

for kk=1:NN
  
  % Unpack arguments
  MU = f.C; s = f.s; W = f.W(:,1:end-1)'; w = f.W(:,end)';
  N = size(X2,1); K = size(X2,2)/2; P = size(Y2,2)/2; M = size(W,1);
  Ax = (Ab(kk,1:4*K))'; bx = (Ab(kk,4*K+1:6*K))'; 
  Cy = (Ab(kk,6*K+1:6*K+4*P))'; dy = (Ab(kk,6*K+4*P+1:end))';
  
  % Transpose matrices, to conform to math notation
  X = X2'; Y = Y2'; MU = MU'; W = W'; w = w';
  
  % Auxiliary variables
  gAx = zeros(1,4*K); gbx = zeros(1,2*K);	% desired error gradients
  gCy =  zeros(1,4*P); gdy = zeros(1,2*P);	% desired error gradients
  f0 = 0;
  
  if (nargout>2)
    H = zeros((6*K+6*P),(6*K+6*P));
  end
  
  for n=1:N
    % Linearly transformed Y and X
    Xn = X(:,n); Xn2 = reshape(Xn,2,K);
    Yn = Y(:,n);
    
    blkx = zeros(2*K, 4*K);
    for k=1:K
      %blkx(2*k-1:2*k,4*(k-1)+1:4*k) = kron(Xn2(:,k)', eye(2));  % original version
      blkx(2*k-1:2*k,4*(k-1)+1:4*k) = [Xn2(1,k) 0 Xn2(2,k) 0;0 Xn2(1,k) 0 Xn2(2,k)]; % fast version
    end
    
    XX = blkx * Ax + bx; % g(X)
    Phi = exp(-sqdist(MU',XX')/(2*s*s));	% Mx1, phi(g(Xn)) is column n
    WPhi = W * diag(sparse(Phi));
    Ytn = W * Phi + w; % f(g(X))
    
    Ytn2 = reshape(Ytn,2,P);
    blkyt = zeros(2*P, 4*P); blkC = zeros(2*P, 2*P);
    for p=1:P
      %blkyt(2*p-1:2*p,4*(p-1)+1:4*p) = kron(Ytn2(:,p)',eye(2));  % original version
      blkyt(2*p-1:2*p,4*(p-1)+1:4*p) = [Ytn2(1,p) 0 Ytn2(2,p) 0;0 Ytn2(1,p) 0 Ytn2(2,p)]; % fast version
      blkC(2*(p-1)+1:2*p,2*(p-1)+1:2*p) = reshape(Cy(4*(p-1)+1: 4*p),2,2);
    end
    
    tmpC = 1/(s*s) * blkC * WPhi * (XX*ones(1,M) - MU)';  % modification
    
    gv_Cy = -blkyt;
    gv_dy = -eye(2*P);
    %gv_Ax = 1/(s*s) * blkC * WPhi * (XX*ones(1,M) - MU)'*blkx;
    %gv_bx = 1/(s*s) * blkC * WPhi * (XX*ones(1,M) - MU)'*eye(2*K);
    gv_Ax = tmpC*blkx;
    gv_bx = tmpC*eye(2*K);
    
    R = Yn - blkC * Ytn - dy;
    f0 = f0 + R'*R; % error value
    
    gAx = gAx + R'*gv_Ax;
    gbx = gbx + R'*gv_bx;
    gCy = gCy + R'*gv_Cy;
    gdy = gdy + R'*gv_dy;
   
    if (nargout>2)
      H = H + [gv_Ax gv_bx gv_Cy gv_dy]'*[gv_Ax gv_bx gv_Cy gv_dy];
    end
    
  end
  
  gAx = 2*gAx; gbx = 2*gbx; gCy = 2*gCy; gdy = 2*gdy; g = [gAx gbx gCy gdy];
  
  if exist('l','var') & ~isempty(l) & (l~=0)
    cx = 0; cy = 0; D = 2;
    gv_Al = zeros(1,4*K); gv_Cl = zeros(1,4*P); 
    for k=1:K
      A = reshape(Ax(4*(k-1)+1:4*k),2,2);
      nA2 = norm(A,'fro')^2; detA = det(A);
      cx = cx + (nA2-D*(detA^2)^(1/D));
      
      iAt = inv(A)';
      gv_Al(4*(k-1)+1: 4*k) = 2*(A(:)' - ((detA^2)^(1/D))*iAt(:)');
    end
    for p=1:P
      A = reshape(Cy(4*(p-1)+1: 4*p),2,2);
      nA2 = norm(A,'fro')^2; detA = det(A);
      cy = cy + (nA2-D*(detA^2)^(1/D));
      
      iAt = inv(A)';
      gv_Cl(4*(p-1)+1:4*p) = 2*(A(:)' - ((detA^2)^(1/D))*iAt(:)');
    end
    
    f0 = f0 + l*(cx + cy);
    g = g + l*([gv_Al zeros(1,2*K) gv_Cl zeros(1,2*P)]);
    
    if (nargout>2)
      H = H + l*[gv_Al zeros(1,2*K) gv_Cl zeros(1,2*P)]'*[gv_Al zeros(1,2*K) gv_Cl zeros(1,2*P)];
    end
  end
  
  ff = [ff; f0];
  gg = [gg; g];
  
  if (nargout>2)
    HH = [HH; 4*H];
  end
      
end
