% [f,fX,E] = rbftrain(X,Y,C[,s,l]) Train RBF to map f(X) = Y
%
% Examples of use (C is a matrix):
% - rbftrain(X,Y,10,[],1e-5): compute 10 BF centres, width and weights.
% - rbftrain(X,Y,10,2,1e-5): compute 10 BF centres and weights; width=2.
% - rbftrain(X,Y,C,[],1e-5): compute BF centres with kmeans initialised at C;
%   compute width and weights.
% - rbftrain(X,Y,C,[1 2],1e-5): like rbftrain(X,Y,C,[],1e-5) but using [1 2]
%   as guess range for the width.
% - rbftrain(X,Y,C,2,1e-5): BF centres=C; width=2; compute weights.
%
% In:
%   X: NxL matrix, N L-dim data points rowwise.
%   Y: NxD matrix, N D-dim data points rowwise.
%   C: if a scalar, the #BFs to use;
%      else, MxL matrix, M L-dim initial values of the RBF centres rowwise.
%   s: RBF width, either a scalar or a range [s1 s2] as initial guess range
%      for the width. Default: calculate (guessing a range).
%      If C is given as a matrix and s as a scalar, only the weights are fit.
%   l: (nonnegative scalar) regularisation parameter. Default: 0.
% Out:
%   f: (struct) the RBF, with fields:
%      type='rbf', centres C (MxL), width s, weights W (Dx(M+1), incl. bias),
%      regularisation parameter l.
%   fX: NxD matrix, f(X).
%   E: 1x2 list, fit error and regularisation error.

% Copyright (c) 2009 by Miguel A. Carreira-Perpinan

function [f,fX,E] = rbftrain(X,Y,C,s,l)

% ---------- Argument defaults ----------
if ~exist('s','var') s = []; end;
if (length(s)==2) & (s(1)==s(2)) s = s(1); end;
if ~exist('l','var') l = []; f.l = 0; else f.l = l; end;
% ---------- End of "argument defaults" ----------

global sqdCX

N = size(X,1); f.type = 'rbf';

% Fit centres
if isscalar(C)
  M = C; e = Inf;	% run kmeans 10 times with rnd init and pick best
  for i=1:10
    [tmpC,tmp,tmpe] = kmeans(X,M);
    if tmpe(end) < e C = tmpC; e = tmpe(end); end
  end
elseif isempty(s)
  C = kmeans(X,size(C,1),C);
end
f.C = C;

% Fit width, if not given
if length(s)~=1
  % Train on 70% of the data and testing on the remaining 30%...
  tmp = randperm(N); tmpT = tmp(1:floor(0.7*N)); tmpt= tmp(length(tmpT)+1:end);
  XT = X(tmpT,:); YT = Y(tmpT,:); Xt = X(tmpt,:); Yt = Y(tmpt,:); Eold = Inf;
  sqdCX = sqdist(C,XT);
  % ...and try log-grid around guess interval for s
  if isempty(s)
    % Guess interval around s=2*mean distance from centre to 4 nearest centres
    tmp = sort(sqdist(C),2); s = 2*mean(mean(sqrt(tmp(:,1:4)),2));
    s = [s/3 3*s]; % s = [s/2 2*s];
  end
  sgrid = exp(linspace(log(s(1)),log(s(2)),17));
  warning('off','MATLAB:nearlySingularMatrix');
  bracketed = 0; decr = 1; incr = 1;
  while ~bracketed
    % Pick width that produces the smallest fit error on the test subset
    for ss = sgrid
      [W,w] = rbftrainW(XT,YT,C,ss,l); f.s = ss; f.W = [W w];
      E = sum(sum((Yt-rbf(Xt,f)).^2));	% Only fit (no regularisation) error
      fprintf('s=%f, E=%f\n',ss,E);
      if E < Eold s = ss; Eold = E; end
    end
    % Enlarge guess interval if the minimum lies on the boundary
    if decr & (s==sgrid(1))
      sgrid = exp(linspace(log(0.99*sgrid(1)/2),log(0.99*sgrid(1)),8));
      incr = 0;
    elseif incr & (s==sgrid(end))
      sgrid = exp(linspace(log(1.01*sgrid(end)),log(2*1.01*sgrid(end)),8));
      decr = 0;
    else bracketed = 1;
    end
  end
  warning('on','MATLAB:nearlySingularMatrix');
end
f.s = s;

% Fit weights
sqdCX = sqdist(C,X);
if nargout>2 [W,w,fX,E] = rbftrainW(X,Y,C,s,l);
elseif nargout>1 [W,w,fX] = rbftrainW(X,Y,C,s,l);
else [W,w] = rbftrainW(X,Y,C,s,l);
end
f.W = [W w];


function [W,w,fX,E] = rbftrainW(X,Y,C,s,l)

global sqdCX

N = size(X,1); M = size(C,1);
G = exp(-sqdCX/(2*s*s)); G1 = sum(G,2); GG = G*G' - G1*(G1'/N);
if ~isempty(l) GG = GG + spdiags(repmat(l,M,1),0,M,M); end
W = (Y'*G'-mean(Y',2)*G1') / GG; w = mean(Y'-W*G,2);

if nargout>2 f.type='rbf'; f.C=C; f.s=s; f.W=[W w]; fX = rbf(X,f); end
if nargout>3
  E = sum(sum((Y-fX).^2));					% Fit error
  if ~isempty(l) E = [E l*(W(:)'*W(:))]; else E = [E 0]; end	% Regul. error
end

