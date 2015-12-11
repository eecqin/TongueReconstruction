% [ff,gg,HH] = Ferr2_linf2(Ab,f,X,Y,l) Local adaptation error with linear function
%
% Reconstruction error E(A,b) wrt the linear transformation (given by A (2x2)
% and b (2x1)), for a linear predictive mapping.
%
% In:
%   Ab: 1x6 vector = [A(:);b(:)]'.
%   f: struct with fields W, w for the linear mapping (see linftrain.m).
%   X: Nx(2.K) matrix containing N 2K-dim contour points rowwise.
%   Y: Nx(2.P) matrix containing N 2P-dim contour points rowwise.
%   l: regularisation parameter to minimise cond(A). Default: 0.
% Out:
%   ff: error value.

% Copyright (c) 2010 by Chao Qin and Mohsen Farhadloo

function [ff,gg,HH] = Ferr2_linf2(Ab,f,X2,Y2,l)

NN = size(Ab,1);
ff = []; gg = []; HH = [];

for kk=1:NN
  
  % Unpack arguments
  W = f.W; w = f.w;
  N = size(X2,1); K = size(X2,2)/2; P = size(Y2,2)/2; 
  Ax = (Ab(kk,1:4*K))'; bx = (Ab(kk,4*K+1:6*K))'; 
  Cy = (Ab(kk,6*K+1:6*K+4*P))'; dy = (Ab(kk,6*K+4*P+1:end))';
  
  % Transpose matrices, to conform to math notation
  X = X2'; Y = Y2';
  
  % Auxiliary variables
  gAx = zeros(1,4*K); gbx = zeros(1,2*K);	% desired error gradients
  gCy =  zeros(1,4*P); gdy = zeros(1,2*P);	% desired error gradients
  f0 = 0;
  H = zeros((6*K+6*P),(6*K+6*P));
  
  for n=1:N
    % Linearly transformed Y and X
    Xn = X(:,n); Xn2 = reshape(Xn,2,K);
    Yn = Y(:,n);
    
    blkx = zeros(2*K, 4*K);
    for k=1:K
      %blkx(2*k-1:2*k,4*(k-1)+1:4*k) = kron(Xn2(:,k)',eye(2)); % original version
      blkx(2*k-1:2*k,4*(k-1)+1:4*k) = [Xn2(1,k) 0 Xn2(2,k) 0;0 Xn2(1,k) 0 Xn2(2,k)]; % fast version
    end
    
    XX = blkx*Ax + bx;   % g(X)
    Ytn = W*XX + w;      % f(g(X))
    
    Ytn2 = reshape(Ytn,2,P);
    blkyt = zeros(2*P, 4*P); blkC = zeros(2*P, 2*P);
    for p=1:P
      %blkyt(2*p-1:2*p,4*(p-1)+1:4*p) = kron(Ytn2(:,p)',eye(2));
      blkyt(2*p-1:2*p,4*(p-1)+1:4*p) = [Ytn2(1,p) 0 Ytn2(2,p) 0;0 Ytn2(1,p) 0 Ytn2(2,p)]; % fast version
      blkC(2*(p-1)+1:2*p,2*(p-1)+1:2*p) = reshape(Cy(4*(p-1)+1:4*p),2,2);
    end
    
    gv_Cy = -blkyt;
    gv_dy = -eye(2*P);
    gv_Ax = -blkC * W * blkx;
    gv_bx = -blkC * W * eye(2*K);
    
    R = Yn - blkC * Ytn - dy;
    f0 = f0 + R'*R;         % error value
    
    gAx = gAx + R'*gv_Ax;
    gbx = gbx + R'*gv_bx;
    gCy = gCy + R'*gv_Cy;
    gdy = gdy + R'*gv_dy;
    
    H = H + [gv_Ax gv_bx gv_Cy gv_dy]'*[gv_Ax gv_bx gv_Cy gv_dy];
    
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
    H = H + l*[gv_Al zeros(1,2*K) gv_Cl zeros(1,2*P)]'*[gv_Al zeros(1,2*K) gv_Cl zeros(1,2*P)];
  end
  
  ff = [ff; f0];
  gg = [gg; g];
  HH = [HH; 4*H];
  
end

