% [ff,gg,HH] = Ferr_rbf(Ab,f,X,Y,l) Global adaptation error with RBF mapping
%
% Reconstruction error E(A,b) wrt the linear transformation (given by A (2x2)
% and b (2x1)), for an RBF predictive mapping:
% - RBF mapping y = f(x) = W.phi(x)+w, phi(x)_m = exp(-|(x-mu_m)/s|²/2)
% - 2D point transformation: xx = g(x) = A.x+b, yy = g(y) = A.y+b
% - Reconstruction error E(A,b) = \sum_n{|g(yn) - f(g(xn))|²}
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

function [ff,gg,HH] = Ferr_rbf(Ab,f,X,Y,l)

% Unpack arguments
A = reshape(Ab(1:4),2,2); b = Ab(5:6)';
MU = f.C; s = f.s; W = f.W(:,1:end-1)'; w = f.W(:,end)';

N = size(X,1); K = size(X,2)/2; P = size(Y,2)/2; M = size(W,1);

% Transpose matrices, to conform to math notation
X = X'; Y = Y'; MU = MU'; W = W'; w = w';

% Linearly transformed Y and X
YY = bsxfun(@plus,kron(eye(P),A)*Y,kron(ones(P,1),b));		% g(Y)
XX = bsxfun(@plus,kron(eye(K),A)*X,kron(ones(K,1),b));		% g(X)

% Auxiliary variables
Phi = exp(-sqdist(MU',XX')/(2*s*s));	% MxN, phi(g(Xn)) is column n
WPhi = W*Phi;				% (2.P)xN
I2A = kron(eye(2),A);
meanM = zeros(2,M); for i=1:K meanM=meanM+MU(2*i-1:2*i,:); end; meanM=meanM/K;
R = bsxfun(@minus,YY-WPhi,w);
ff = R(:)'*R(:);			% error value
if exist('l','var') & ~isempty(l)
  nA2 = norm(A,'fro')^2; detA = det(A); D = 2;
  ff = ff + l*(nA2-D*(detA^2)^(1/D));
end

if nargout>1

  gA = zeros(1,4); gb = zeros(1,2);	% desired error gradients
  HH = zeros(6,6);			% desired error approximate Hessian
  
  for n=1:N
    
    Xn = X(:,n); Xn2 = reshape(Xn,2,K);
    Yn = Y(:,n); Yn2 = reshape(Yn,2,P);
    
    meanx = mean(Xn2,2);
    meanKRONxx = mean([Xn2(1,:).^2;prod(Xn2,1);prod(Xn2,1);Xn2(2,:).^2],2);
    
    WPhiMt = (W*(diag(sparse(Phi(:,n)))*MU'))';
    t=zeros(4,2*P);
    for i=1:K t=t+kron(Xn2(:,i),WPhiMt(2*i-1:2*i,:)); end; t=t/K;
    
    gv_A = kron(Yn2',eye(2)) + ...
           K/(s*s)*(WPhi(:,n)*(I2A*meanKRONxx+kron(meanx,b))' - t');
    
    gv_b = kron(ones(P,1),eye(2)) + ...
           K/(s*s)*(WPhi(:,n)*(A*meanx+b)'-W*(diag(sparse(Phi(:,n)))*meanM'));
    
    gA = gA + R(:,n)'*gv_A;
    gb = gb + R(:,n)'*gv_b;
    HH = HH + [gv_A gv_b]'*[gv_A gv_b];
    
  end
  
  if exist('l','var') & ~isempty(l)
%    gv_Al = l*(A(:)' - sign(detA)*[A(2,2) -A(1,2) -A(2,1) A(1,1)]); % D=2
    iAt = inv(A)'; gv_Al = l*(A(:)' - ((detA^2)^(1/D))*iAt(:)');
    gA = gA + gv_Al;
    HH = HH + [gv_Al 0 0]'*[gv_Al 0 0];
  end
  
  gA = 2*gA; gb = 2*gb; gg = [gA gb];
  HH = 4*HH;

end

