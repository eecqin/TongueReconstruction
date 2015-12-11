function [Y,fgX,gX] = adapted_f2(X,f,Ab)  % for fcn F

N = size(X,1); K = size(X,2)/2; P = (length(Ab)-6*K)/6;
fcn = str2func(f.type);			% Use string as function handle

% Unpack arguments
Ax = (Ab(1:4*K))'; bx = (Ab(4*K+1:6*K))'; 
Ay = (Ab(6*K+1:6*K+4*P))'; by = (Ab(6*K+4*P+1:end))';

Y = zeros(N,2*P); fgX = Y; gX = zeros(N,2*K);
for n=1:N
  % Linearly transformed Y and X
  Xn = X(n,:)'; Xn2 = reshape(Xn,2,K);
  blkx = zeros(2*K,4*K);
  for k=1:K
    %blkx(2*k-1:2*k,4*(k-1)+1:4*k) = kron(Xn2(:,k)', eye(2));  % original version
    blkx(2*k-1:2*k,4*(k-1)+1:4*k) = [Xn2(1,k) 0 Xn2(2,k) 0;0 Xn2(1,k) 0 Xn2(2,k)]; % fast version
  end
  gX(n,:) = (blkx*Ax + bx)'; % g(X(n))  
  fgX(n,:) = fcn(gX(n,:),f); % f(g(X(n)))
  
  blkiA = zeros(2*P, 2*P);
  for p=1:P
    blkiA(2*(p-1)+1:2*p,2*(p-1)+1:2*p) = inv(reshape(Ay(4*(p-1)+1: 4*p),2,2));
  end
  
  Y(n,:) = (blkiA * fgX(n,:)' - blkiA * by)';
end
