function [ff,gg,HH] = Ferr2_rbf2a(Ab,f,X2,Y2,Cd)
% $$$ function [ff,gg,HH] = Ferr2_rbf2a(Ab,f,X2,Y2) % debug: for validating gradients only
  
NN = size(Ab,1);
ff = []; gg = []; HH = [];

for kk=1:NN
  
  % Unpack arguments
  MU = f.C; s = f.s; W = f.W(:,1:end-1)'; w = f.W(:,end)';
  N = size(X2,1); K = size(X2,2)/2; P = size(Y2,2)/2; M = size(W,1);
  Ax = (Ab(kk,1:4*K))'; bx = (Ab(kk,4*K+1:6*K))'; 
  Cy = (Cd(kk,1:4*P))'; dy = (Cd(kk,4*P+1:6*P))';
  % $$$ Cy = repmat([1 0 0 1]',P,1); dy = repmat([0 0]',P,1); % debug: for validating gradient only
  
  % Transpose matrices, to conform to math notation
  X = X2'; Y = Y2'; MU = MU'; W = W'; w = w';
  
  % Auxiliary variables
  gAx = zeros(1,4*K); gbx = zeros(1,2*K);	% desired error gradients
  f0 = 0;
  H = zeros(6*K,6*K);
  
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
    
    blkC = zeros(2*P, 2*P);
    for p=1:P
      blkC(2*(p-1)+1:2*p,2*(p-1)+1:2*p) = reshape(Cy(4*(p-1)+1: 4*p),2,2);
    end

    gv_Ax = 1/(s*s) * blkC * WPhi * (XX*ones(1,M) - MU)'*blkx;
    gv_bx = 1/(s*s) * blkC * WPhi * (XX*ones(1,M) - MU)'*eye(2*K);
    
    R = Yn - blkC * Ytn - dy;
    f0 = f0 + R'*R; % error value
    
    gAx = gAx + R'*gv_Ax;
    gbx = gbx + R'*gv_bx;
    
    H = H + [gv_Ax gv_bx]'*[gv_Ax gv_bx];
    
  end
  
  gAx = 2*gAx; gbx = 2*gbx; g = [gAx gbx];
  
  ff = [ff; f0];
  gg = [gg; g];
  HH = [HH; 4*H];
  
end
