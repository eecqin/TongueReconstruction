% [gA,gb] = Fgrad(X,Y,W,v,MU,s,A,b) Adaptation (RBF prediction): error gradient
%
% Gradient of the reconstruction error F wrt the linear transformation (given
% by A (2x2) and b (2x1)), for a RBF prediction:
% - RBF mapping y = f(x) = W.phi(x)+v, phi(x)_m = exp(-|(x-mu_m)/s|²/2)
% - 2D point transformation: xx = g(x) = A.x+b, yy = g(y) = A.y+b
% - Reconstruction error F(A,b) = \sum_n{|g(yn) - f(g(xn))|²}
%
% In:
%   X: Nx(2.K) matrix containing N 2K-dim contour points rowwise.
%   Y: Nx(2.P) matrix containing N 2P-dim contour points rowwise.
%   W: Mx(2.P) matrix containing M 2P-dim weight vectors rowwise.
%   v: 1x(2.P) vector containing the 2P-dim bias vector rowwise.
%   MU: Mx(2.K) matrix containing M 2K-dim BF centres rowwise.
%   s: scalar containing the BF width.
%   A, b: 2x2 matrix and 2x1 vector at which to evaluate the gradient.
% Out:
%   gA: 1x4 gradient wrt A(:)
%   gb: 1x2 gradient wrt b(:)

% Copyright (c) 2008 by Miguel A. Carreira-Perpinan

function [gA,gb] = Fgrad(X,Y,W,v,MU,s,A,b)

N = size(X,1); K = size(X,2)/2; P = size(Y,2)/2; M = size(W,1);

% Transpose matrices, to conform to math notation
X = X'; Y = Y'; MU = MU'; W = W'; v = v';

% Linearly transformed Y and X
YY = kron(eye(P),A)*Y + repmat(kron(ones(P,1),b),1,N);		% g(Y)
XX = kron(eye(K),A)*X + repmat(kron(ones(K,1),b),1,N);		% g(X)

% Auxiliary variables
Phi = exp(-sqdist(MU',XX')/(2*s*s));	% MxN, phi(g(Xn)) is column n
WPhi = W*Phi;				% (2.P)xN
I2A = kron(eye(2),A);
meanM = zeros(2,M); for i=1:K meanM=meanM+MU(2*i-1:2*i,:); end; meanM=meanM/K;
R = YY - WPhi - repmat(v,1,N);

gA = zeros(1,4); gb = zeros(1,2);	% Desired gradients

for n=1:N
  
  Xn = X(:,n); Xn2 = reshape(Xn,2,K);
  Yn = Y(:,n); Yn2 = reshape(Yn,2,P);

  meanx = mean(reshape(Xn,2,K),2);
  meanKRONxx = mean([Xn2(1,:).^2;prod(Xn2,1);prod(Xn2,1);Xn2(2,:).^2],2);

  WPhiMt = (W*diag(Phi(:,n))*MU')';
  t=zeros(4,2*P); for i=1:K t=t+kron(Xn2(:,i),WPhiMt(2*i-1:2*i,:)); end; t=t/K;
  
  gv_A = kron(Yn2',eye(2)) + K/(s*s)*(WPhi(:,n)*(I2A*meanKRONxx+kron(meanx,b))' - t');
  
  gv_b = kron(ones(P,1),eye(2)) + K/(s*s)*(WPhi(:,n)*(A*meanx+b)' - W*diag(Phi(:,n))*meanM');
  
  gA = gA + R(:,n)'*gv_A;
  gb = gb + R(:,n)'*gv_b;
  
end

gA = 2*gA; gb = 2*gb;

