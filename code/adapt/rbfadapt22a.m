function [AbCd,E] = rbfadapt22a(X,Y,f,AbCd,tol,maxit) % for fcn E using alternative opt

% ---------- Argument defaults ----------
% $$$ if ~exist('A','var') | isempty(A) | ~exist('b','var') | isempty(b)
% $$$   A = eye(2); b = [0;0];	% Initial A, b
% $$$ end;
% $$$ if ~exist('l','var') l = []; end;
% $$$ if ~exist('tol','var') | isempty(tol) tol = 1e-5; end;
% ---------- End of "argument defaults" ----------

N = size(X,1); K = size(X,2)/2; P = size(Y,2)/2;
it = 0;
AbCd = AbCd';
switch f.type
 case 'rbf'
  [ff,g,H] = Ferr2_rbf2(AbCd,f,X,Y,0); 
 case 'linf'
  [ff,g,H] = Ferr2_linf2(AbCd,f,X,Y,0); 
end
F = ff;
ff_old = 1e+5;

fprintf(1,'ff=%f  norm(g)=%f\n',ff,norm(g));

cnt = 0; max_cnt = 5;

while (norm(g) > tol) & (it < maxit)
%while (abs(ff-ff_old) > tol*abs(ff_old)) & (it < maxit)
   
  ff_old = ff;
  g_old = g; 
  
  % Unpack arguments
  Ax = AbCd(1:4*K); bx = AbCd(4*K+1:6*K);
  Cy = AbCd(6*K+1:6*K+4*P); dy = AbCd(6*K+4*P+1:end);
  
  gX = zeros(N,2*K); fgX = zeros(N,2*P);
  
% $$$   % 1. Fix (Ax,bx) obtained from rbfadapt.m and solve for (Cy,dy) by least square
% $$$   [Cy,dy] = linfadapt22a(X,Y,f,Ax,bx);
% $$$   Ab = [Ax bx]; Cd = [Cy dy];
% $$$   
% $$$   % 2. Fix (Cy,dy) and solve for (Ax,bx) by BFGS using Ferr2_rbf2a
% $$$   switch f.type
% $$$    case 'rbf'
% $$$     [tmp1,tmp2,H0] = Ferr2_rbf2a(Ab,f,X,Y,Cd); % Initial Hessian estimate
% $$$     [tmp1,tmp2] = Obfgs(@Ferr2_rbf2a,{f,X,Y,Cd},Ab',tol,[],[],inv(H0));
% $$$    case 'linf'
% $$$     [tmp1,tmp2,H0] = Ferr2_linf2a(Ab,f,X,Y,Cd); % Initial Hessian estimate
% $$$     [tmp1,tmp2] = Obfgs(@Ferr2_linf2a,{f,X,Y,Cd},Ab',tol,[],[],inv(H0));
% $$$   end
% $$$   Ax = tmp1(end,1:4*K); bx = tmp1(end,4*K+1:6*K);
  
  % 1. Fix (Cy,dy) obtained from rbfadapt and solve for (Ax,bx) by BFGS using Ferr2_rbf2a
  Ab = [Ax bx]; Cd = [Cy dy];
  switch f.type
   case 'rbf'
    [tmp1,tmp2,H0] = Ferr2_rbf2a(Ab,f,X,Y,Cd); % Initial Hessian estimate
    [tmp1,tmp2] = Obfgs(@Ferr2_rbf2a,{f,X,Y,Cd},Ab',tol,[],[]);
   case 'linf'
    [tmp1,tmp2,H0] = Ferr2_linf2a(Ab,f,X,Y,Cd); % Initial Hessian estimate
    [tmp1,tmp2] = Obfgs(@Ferr2_linf2a,{f,X,Y,Cd},Ab',tol,[],[]);
  end
  Ax = tmp1(end,1:4*K); bx = tmp1(end,4*K+1:6*K);
  
  % 2. Fix (Ax,bx) and solve for (Cy,dy) by least square using linftrain.m
  [Cy,dy] = linfadapt22a(X,Y,f,Ax,bx);
  
  
  AbCd = [Ax bx Cy dy];
  switch f.type
   case 'rbf'
    [ff,g,H] = Ferr2_rbf2(AbCd,f,X,Y,0);
   case 'linf'
    [ff,g,H] = Ferr2_linf2(AbCd,f,X,Y,0);
  end
  
  %if norm(g)==norm(g_old)  % original
  if (abs(norm(g)-norm(g_old))<1e-5)
    cnt = cnt + 1;
    if (cnt>max_cnt) break; end;
  else
    cnt = 0;    % reset
  end
  
  F = [F; ff];	% Update lists
  it = it + 1;
  
  fprintf(1,'ff=%f  norm(g)=%f\n',ff,norm(g));
end

E = F(end);

