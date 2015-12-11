% C = pcorrel(A,B)
% Compute Pearson's correlation between rowwise data sets A and B
%
% In:
%  See rmse.
% Out:
%  C: Dx1 vector of Pearson's correlation for each articulatory channel.

% Copyright (c) 2008 by Chao Qin and Miguel A. Carreira-Perpinan

function C = pcorrel(A,B)
  
[N,D] = size(A); 

C = zeros(D,1);

% Centre the dataset
A = A - repmat(mean(A,1),N,1); B = B - repmat(mean(B,1),N,1);

% Vectorised version to compute the Pearson correlation
C = (sum(A.*B,1)./(sqrt(sum(A.^2,1)).*sqrt(sum(B.^2,1))))';

