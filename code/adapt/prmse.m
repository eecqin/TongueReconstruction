% [E,Estd] = prmse(A,B[,scale]) Pointwise root-mean-square error (L2-norm, ie. Euclidean norm) 
% between estimation and measurements of tongue contours

% In:
%  A: NxD matricies of estimations
%  B: NxD matricies of ground truth measurements
%  scale: scale between pixel and mm
% 
% Out:
%  E: D/2x1 vector containing D/2 pointwise errors
%  Estd: D/2x1 vector containing D/2 pointwise error standard deviations

% Copyright (c) 2010 by Chao Qin

function [E,Estd] = prmse(A,B,scale)

% ------ Argument default --------------------------------
if ~exist('scale','var') | isempty(scale) scale = 1; end;
% ------ End of argument default -------------------------

[N,D] = size(A); 

% Scaling
A = A/scale; B = B/scale;

E = zeros(D/2,1); Estd = E;
for i=1:D/2
  E(i) = sqrt(mean( sqrt((A(:,2*i-1)-B(:,2*i-1)).^2+(A(:,2*i)-B(:,2*i)).^2) ));
  Estd(i) = sqrt(std( sqrt((A(:,2*i-1)-B(:,2*i-1)).^2+(A(:,2*i)-B(:,2*i)).^2) ));
end

