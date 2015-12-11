% [ff,gg,HH] = Ferr_rbf2(Ab,f,X,Y,l) Adaptation error with RBF mapping
%
% Reconstruction error E(A,b) wrt the linear transformation (given by A (2x2)
% and b (2x1)), for an RBF predictive mapping:
% - RBF mapping y = f(x) = W.phi(x)+w, phi(x)_m = exp(-|(x-mu_m)/s|�/2)
% - 2D point transformation: xx = g(x) = A.x+b, yy = g(y) = A.y+b
% - Reconstruction error E(A,b) = \sum_n{|g(yn) - f(g(xn))|�}
%
% In:
%   Ab: 1x6 vector = [A(:);b(:)]'.
%   f: struct with fields C, s, W for the RBF (see rbftrain.m).
%   X: Nx(2.K) matrix containing N 2K-dim contour points rowwise.
%   Y: Nx(2.P) matrix containing N 2P-dim contour points rowwise.
%   l: regularisation parameter to minimise cond(A). Default: 0.
% Out:
%   ff: error value.
%   gg: 1x6 error gradient.
%   HH: 6x6 error approximate Hessian (a la Gauss-Newton).

% Copyright (c) 2009 by Miguel A. Carreira-Perpinan

% K = 3; P = 20; N = 100; M = 20;
% x = linspace(0,2,P); y = sin((1+randn(N,1)/5)*x+repmat(randn(N,1)/10,1,P));
% x = repmat(x,N,1)+randn(N,P)/40;
% Y1 = [];for i=1:N tmp = [x(i,:);y(i,:)]; Y1 = [Y1;tmp(:)']; end; Y1 = 70*Y1;
% X1 = Y1(:,[5 6 11 12 31 32]); % K = 3;
%
% % 2D-wise transformation
% iA0 = 0.5*[cos(0.7) -sin(0.7);sin(0.7) cos(0.7)]; ib0 = [-100;50];
% A0 = inv(iA0); b0 = -inv(iA0)*ib0; % This is the ground truth
% Y2 = (kron(eye(P),iA0)*Y1')'+repmat(repmat(ib0',1,P),N,1);
% X2 = Y2(:,[5 6 11 12 31 32]);
%
% l = 0; [f,fX,E] = rbftrain(X1,Y1,M,[],l);	% Fit and prediction error
% tol = 1e-1;					% Set this carefully
%
% A = eye(2); b = [0;0];	% Initial A, b
% Ab = repmat ([A(:); b(:)],K+P,1);

function [ff,gg,HH] = Ferr_rbf2(Ab,f,X2,Y2,l)

NN = size(Ab,1);
ff = []; gg = []; HH = [];

for kk=1:NN
  
  % Unpack arguments
  MU = f.C; s = f.s; W = f.W(:,1:end-1)'; w = f.W(:,end)';
  N = size(X2,1); K = size(X2,2)/2; P = size(Y2,2)/2; M = size(W,1);
  Ax = (Ab(kk,1:4*K))'; bx = (Ab(kk,4*K+1:6*K))'; 
  Ay = (Ab(kk,6*K+1:6*K+4*P))'; by = (Ab(kk,6*K+4*P+1:end))';
  
  % Transpose matrices, to conform to math notation
  X = X2'; Y = Y2'; MU = MU'; W = W'; w = w';
  
  % Auxiliary variables
  gAx = zeros(1,4*K); gbx = zeros(1,2*K);	% desired error gradients
  gAy =  zeros(1,4*P); gby = zeros(1,2*P);	% desired error gradients  
  f0 = 0;
  H = zeros((6*K+6*P),(6*K+6*P));
  
  for n=1:N
    % Linearly transformed Y and X
    Xn = X(:,n); Xn2 = reshape(Xn,2,K);
    Yn = Y(:,n); Yn2 = reshape(Yn,2,P);
    blkx = zeros(2*K, 4*K);
    blky = zeros(2*P, 4*P);
    for k=1:K
      %blkx(2*k-1:2*k,4*(k-1)+1:4*k) = kron(Xn2(:,k)', eye(2));  % original version
      blkx(2*k-1:2*k,4*(k-1)+1:4*k) = [Xn2(1,k) 0 Xn2(2,k) 0;0 Xn2(1,k) 0 Xn2(2,k)]; % fast version
    end
    
    for p=1:P
      %blky(2*p-1:2*p,4*(p-1)+1:4*p) = kron(Yn2(:,p)', eye(2));  % original version
      blky(2*p-1:2*p,4*(p-1)+1:4*p) = [Yn2(1,p) 0 Yn2(2,p) 0;0 Yn2(1,p) 0 Yn2(2,p)]; % fast version
    end
    
    XX = blkx * Ax + bx; % g(X)
    YY = blky * Ay + by; % g(Y)
    
    Phi = exp(-sqdist(MU',XX')/(2*s*s));	% Mx1, phi(g(Xn)) is column n
    WPhi = W*diag(sparse(Phi));			
    gv_Ax = 1/(s*s) * WPhi * (XX*ones(1,M) - MU)'*blkx;
    gv_bx = 1/(s*s) * WPhi * (XX*ones(1,M) - MU)'*eye(2*K);
    gv_Ay = blky;
    gv_by = eye(2*P);
    
    R = YY - W * Phi - w;
    f0 = f0 + R'*R;               % error value
    
    gAx = gAx + R'*gv_Ax;
    gbx = gbx + R'*gv_bx;
    gAy = gAy + R'*gv_Ay;
    gby = gby + R'*gv_by;
    
    H = H + [gv_Ax gv_bx gv_Ay gv_by]'*[gv_Ax gv_bx gv_Ay gv_by];
    
  end
  
  gAx = 2*gAx; gbx = 2*gbx; gAy = 2*gAy; gby = 2*gby; g = [gAx gbx gAy gby];
  ff = [ff; f0];
  gg = [gg; g];
  HH = [HH; 4*H];
  
end
