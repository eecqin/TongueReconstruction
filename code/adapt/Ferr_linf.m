% ff = Ferr_linf(Ab,f,X,Y,l) Global adaptation error with linear function
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

% Copyright (c) 2009 by Miguel A. Carreira-Perpinan

function ff = Ferr_linf(Ab,f,X,Y,l)

% Unpack arguments
A = reshape(Ab(1:4),2,2); b = Ab(5:6)';
W = f.W; w = f.w;

K = size(X,2)/2; P = size(Y,2)/2;

% Transpose matrices, to conform to math notation
X = X'; Y = Y';

% Linearly transformed Y and X
YY = bsxfun(@plus,kron(eye(P),A)*Y,kron(ones(P,1),b));		% g(Y)
XX = bsxfun(@plus,kron(eye(K),A)*X,kron(ones(K,1),b));		% g(X)

R = bsxfun(@minus,YY-W*XX,w);
ff = R(:)'*R(:);
if exist('l','var') & ~isempty(l)
  nA2 = norm(A,'fro')^2; detA = det(A); D = 2;
  ff = ff + l*(nA2-D*(detA^2)^(1/D));
end

% $$$ % Another method
% $$$ R = zeros(2*P,N);
% $$$ Q = kron(ones(P,1),eye(2)); for i=1:K Q = Q - W(:,2*i-1:2*i); end;
% $$$ for n=1:N
% $$$   Xn = X(:,n); Xn2 = reshape(Xn,2,K);
% $$$   Yn = Y(:,n); Yn2 = reshape(Yn,2,P);
% $$$   PPn = kron(Yn2',eye(2));
% $$$   for i=1:K PPn = PPn - kron(Xn2(:,i)',W(:,2*i-1:2*i)); end;
% $$$   R(:,n) = PPn*A(:) + Q*b(:) - w;
% $$$ end
% $$$ ff = R(:)'*R(:);

