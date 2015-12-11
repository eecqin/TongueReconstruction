function [ff,gg] = Ferr_linf2(Ab,f,X2,Y2)
  
NN = size(Ab,1);
ff = []; gg = [];

for kk=1:NN
  
  % Unpack arguments
  W = f.W; w = f.w;
  N = size(X2,1); K = size(X2,2)/2; P = size(Y2,2)/2; M = size(W,1);  
  Ax = (Ab(kk,1:4*K))'; bx = (Ab(kk,4*K+1:6*K))'; 
  Ay = (Ab(kk,6*K+1:6*K+4*P))'; by = (Ab(kk,6*K+4*P+1:end))';
  
  % Transpose matrices, to conform to math notation
  X = X2'; Y = Y2'; 
  
  % Auxiliary variables
  gAx = zeros(1,4*K); gbx = zeros(1,2*K);	% desired error gradients
  gAy =  zeros(1,4*P); gby = zeros(1,2*P);	% desired error gradients
  f0 = 0;
  
  for n=1:N
    % Linearly transformed Y and X
    Xn = X(:,n); Xn2 = reshape(Xn,2,K);
    Yn = Y(:,n); Yn2 = reshape(Yn,2,P);
    blkx = zeros(2*K, 4*K);
    blky = zeros(2*P, 4*P);
    for k=1:K
      %blkx(2*k-1:2*k,4*(k-1)+1:4*k) = kron(Xn2(:,k)', eye(2));  % original version
      blkx(2*k-1:2*k,4*(k-1)+1:4*k) = [Xn2(1,k) 0 Xn2(2,k) 0;0 Xn2(1,k) 0 Xn2(2,k)];  % fast version
    end
    
    for p=1:P
      %blky(2*p-1:2*p,4*(p-1)+1:4*p) = kron(Yn2(:,p)', eye(2));  % original version
      blky(2*p-1:2*p,4*(p-1)+1:4*p) = [Yn2(1,p) 0 Yn2(2,p) 0;0 Yn2(1,p) 0 Yn2(2,p)];  % fast version
    end
    
    XX = blkx * Ax + bx;      % g(X)
    YY = blky * Ay + by;      % g(Y)
    
    Pnx = -W * blkx;
    Qx = -W;
    Pny = blky;
    Qy = eye(2*P);
    
    gv_Ax = Pnx;
    gv_bx = Qx;
    gv_Ay = Pny;
    gv_by = Qy;
    
    R = YY - W * XX - w;
    f0 = f0 + R'*R;          % error value
    
    gAx = gAx + R'*gv_Ax;
    gbx = gbx + R'*gv_bx;
    gAy = gAy + R'*gv_Ay;
    gby = gby + R'*gv_by;
    
  end
  
  gAx = 2*gAx; gbx = 2*gbx; gAy = 2*gAy; gby = 2*gby; g = [gAx gbx gAy gby];
  ff = [ff; f0];
  gg = [gg; g];
  
end