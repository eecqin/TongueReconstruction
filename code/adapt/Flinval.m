% F = Flinval(X,Y,W,v,A,b) Adaptation (linear prediction): error value
%
% Value of the reconstruction error F wrt the linear transformation (given
% by A (2x2) and b (2x1)), for a linear prediction mapping: see Flin.
%
% In: see Flin.
% Out:
%   F: scalar containing the value of F(A,b).

% Copyright (c) 2008 by Miguel A. Carreira-Perpinan

function F = Flinval(X,Y,W,v,A,b)

N = size(X,1); K = size(X,2)/2; P = size(Y,2)/2;

% Transpose matrices, to conform to math notation
X = X'; Y = Y'; W = W'; v = v';

% Linearly transformed Y and X
YY = kron(eye(P),A)*Y + repmat(kron(ones(P,1),b),1,N);		% g(Y)
XX = kron(eye(K),A)*X + repmat(kron(ones(K,1),b),1,N);		% g(X)

R = YY - W*XX - repmat(v,1,N);
F = R(:)'*R(:);

% $$$ % Another method
% $$$ R = zeros(2*P,N);
% $$$ Q = kron(ones(P,1),eye(2)); for i=1:K Q = Q - W(:,2*i-1:2*i); end;
% $$$ for n=1:N
% $$$   Xn = X(:,n); Xn2 = reshape(Xn,2,K);
% $$$   Yn = Y(:,n); Yn2 = reshape(Yn,2,P);
% $$$   PPn = kron(Yn2',eye(2));
% $$$   for i=1:K PPn = PPn - kron(Xn2(:,i)',W(:,2*i-1:2*i)); end;
% $$$   R(:,n) = PPn*A(:) + Q*b(:) - v;
% $$$ end
% $$$ F = R(:)'*R(:);

