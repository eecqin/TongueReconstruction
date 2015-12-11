function [Ab,E] = linfadapt22(X,Y,f,Ab,l,tol) % for fcn E

% ---------- Argument defaults ----------
% $$$ if ~exist('A','var') | isempty(A) | ~exist('b','var') | isempty(b)
% $$$   A = eye(2); b = [0;0];	% Initial A, b
% $$$ end;
% $$$ if ~exist('l','var') l = []; end;
% $$$ if ~exist('tol','var') | isempty(tol) tol = 1e-5; end;
% ---------- End of "argument defaults" ----------

  
[tmp1,tmp2,H0] = Ferr2_linf2(Ab',f,X,Y,l); % Initial Hessian estimate
[tmp1,tmp2] = Obfgs(@Ferr2_linf2,{f,X,Y,l},Ab,tol,[],[]);
Ab = tmp1(end,:);
E = tmp2(end);

