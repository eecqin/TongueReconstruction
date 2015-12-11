% [E,Estd] = tonguermse(A,B[,scale]) Root-mean-square error (L2-norm, ie. Euclidean norm) 
% between estimation and measurements of tongue contours

% In:
%  A: NxD matricies of estimations
%  B: NxD matricies of ground truth measurements
%  scale: scale between pixel and mm
% 
% Out:
%  E: Dx1 vector containing D/2 pointwise errors
%  Estd: Dx1 vector containing D/2 pointwise error standard deviations

% Copyright (c) 2007 by Chao Qin

function [E,Estd] = tonguermse(A,B,scale)

% ------ Argument default --------------------------------
if ~exist('scale','var') | isempty(scale) scale = 1; end;
% ------ End of argument default -------------------------

[N,D] = size(A); 

% Scaling
A = A/scale; B = B/scale;

% original version
%E = zeros(D,1); Estd = zeros(D,1);
%for i=1:D
%  E(i) = sqrt(mean((A(:,i)-B(:,i)).^2));
%  Estd(i) = sqrt(std((A(:,i)-B(:,i)).^2));
%end

E = zeros(D/2,1); Estd = E;
for i=1:D/2
  E(i) = sqrt(mean( sqrt((A(:,2*i-1)-B(:,2*i-1)).^2+(A(:,2*i)-B(:,2*i)).^2) ));
  Estd(i) = sqrt(std( sqrt((A(:,2*i-1)-B(:,2*i-1)).^2+(A(:,2*i)-B(:,2*i)).^2) ));
end

%E = zeros(D/2,1); 
%for i=1:D/2
%  for n=1:N
%    E(i) = E(i) + norm(A(n,[2*i-1:2*i])-B(n,[2*i-1:2*i]));
%  end
%  E(i) = sqrt(E(i)/N);  
%end
