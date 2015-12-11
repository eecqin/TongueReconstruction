% [f,fX,E] = rbftrain(X,Y,C[,s,l]) Train RBF to map f(X) = Y
%
% Examples of use (C is a matrix):
% - rbftrain(X,Y,10,[],1e-5): compute 10 BF centres, width and weights.
% - rbftrain(X,Y,10,2,1e-5): compute 10 BF centres and weights; width=2.
% - rbftrain(X,Y,C,[],1e-5): compute BF centres with kmeans initialised at C;
%   compute width and weights.
% - rbftrain(X,Y,C,2,1e-5): BF centres=C; width=2; compute weights.
%
% In:
%   X: NxL matrix, N L-dim data points rowwise.
%   Y: NxD matrix, N D-dim data points rowwise.
%   C: if a scalar, the #BFs to use;
%      else, MxL matrix, M L-dim initial values of the RBF centres rowwise.
%   s: scalar RBF width. Default: calculate.
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
if ~exist('l','var') l = []; f.l = 0; else f.l = l; end;
% ---------- End of "argument defaults" ----------

global sqdC sqdCX

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
f.C = C; sqdC = sqdist(C);

% Fit width, if not given
if isempty(s)
  % Try log-grid around s = 2*mean distance from centre to 4 nearest centres...
  tmp = sort(sqdC,2); s = 2*mean(mean(sqrt(tmp(:,1:4)),2));
  sgrid = exp(linspace(log(s/3),log(3*s),21));
%  sgrid = exp(linspace(log(s/2),log(2*s),11));
  % ...training on 70% of the data and testing on the remaining 30%
  tmp = randperm(N); tmpT = tmp(1:floor(0.7*N)); tmpt= tmp(length(tmpT)+1:end);
  XT = X(tmpT,:); YT = Y(tmpT,:); Xt = X(tmpt,:); Yt = Y(tmpt,:); Eold = Inf;
  sqdCX = sqdist(C,XT);
  warning('off','MATLAB:nearlySingularMatrix');
  for ss = sgrid	% Pick width that produces the smallest fit error
    [W,w] = rbftrainW(XT,YT,C,ss,l); f.s = ss; f.W = [W w];
    E = sum(sum((Yt-rbf(Xt,f)).^2));
    fprintf('s=%f, E=%f\n',ss,E);
    if E < Eold s = ss; Eold = E; end
  end
  warning('on','MATLAB:nearlySingularMatrix');
end
f.s = s;

% Fit weights
sqdCX = sqdist(C,X);
if nargout>2
  [W,w,fX,E] = rbftrainW(X,Y,C,s,l);
elseif nargout>1
  [W,w,fX] = rbftrainW(X,Y,C,s,l);
else
  [W,w] = rbftrainW(X,Y,C,s,l);
end
f.W = [W w];


function [W,w,fX,E] = rbftrainW(X,Y,C,s,l)

global sqdC sqdCX

N = size(X,1);
G = exp(-sqdCX/(2*s*s)); G1 = sum(G,2); GG = G*G' - G1*(G1'/N);
if ~isempty(l) lGxx = l*exp(-sqdC/(2*s*s)); GG = GG + lGxx; end
W = (Y'*G'-sum(Y',2)*G1'/N) / GG; w = sum(Y'-W*G,2)/N;

if nargout>2 f.type='rbf'; f.C=C; f.s=s; f.W=[W w]; fX = rbf(X,f); end
if nargout>3
  E = sum(sum((Y-fX).^2));					% Fit error
  if ~isempty(l) E = [E trace(W*lGxx*W')]; else E = [E 0]; end	% Regul. error
end

