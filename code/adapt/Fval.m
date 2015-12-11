% F = Fval(X,Y,W,v,MU,s,A,b) Adaptation (RBF prediction): error value
%
% Value of the reconstruction error F wrt the linear transformation (given
% by A (2x2) and b (2x1)), for a RBF prediction: see gradF.
%
% In: see gradF.
% Out:
%   F: scalar containing the value of F(A,b).

% Copyright (c) 2008 by Miguel A. Carreira-Perpinan

function F = Fval(X,Y,W,v,MU,s,A,b)

N = size(X,1); K = size(X,2)/2; P = size(Y,2)/2; M = size(W,1);

% Transpose matrices, to conform to math notation
X = X'; Y = Y'; MU = MU'; W = W'; v = v';

% Linearly transformed Y and X
YY = kron(eye(P),A)*Y + repmat(kron(ones(P,1),b),1,N);		% g(Y)
XX = kron(eye(K),A)*X + repmat(kron(ones(K,1),b),1,N);		% g(X)

% Auxiliary variables
Phi = exp(-sqdist(MU',XX')/(2*s*s));	% MxN, phi(g(Xn)) is column n
WPhi = W*Phi;				% (2.P)xN
R = YY - WPhi - repmat(v,1,N);
F = R(:)'*R(:);

