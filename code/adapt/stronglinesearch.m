function alpha = stronglinesearch(f,paramf,x,p,ff,g,alpha0,rho,c1,c2)
% Search for a step length to satisfy strong Wolfe condition. 
% Algorithm 3.2 in textbook.

  
% ---------- Argument defaults -----------
if ~exist('ff','var') | isempty(ff) [ff,g] = f(x',paramf{:}); end;
if ~exist('g','var') | isempty(g) [ff,g] = f(x',paramf{:}); end;
if ~exist('alpha0','var') | isempty(alpha0) alpha0 = 1.0; end;
if ~exist('rho','var') | isempty(rho) rho = 0.8; end;
if ~exist('c1','var') | isempty(c1) c1 = 1e-4; end;
if ~exist('c2','var') | isempty(c1) c2 = 0.9; end;  % This value set the quality of line search.
if (g*p>=0) error('linesearch: the search direction is not descent!\n'); end;
% ---------- End of "argument defaults" ----------


amax = 50*alpha0;  % Maximum step length.
maxit = 100;      % Maximum number of iterations.

a0 = 0;
ff0 = ff;
a1 = alpha0;

i = 1;
while 1
  x1 = x + a1*p;
  [ff1,g1] = f(x1',paramf{:});
  if (ff1>ff+c1*a1*g*p) | (ff1>=ff0 & i>1)
    alpha = ZOOM(f,paramf,x,p,ff,g,a0,a1,c1,c2);
    return;
  end
  
  if abs(g1*p)<=-c2*(g*p)
    alpha = a1;
    return;
  end
  
  if (g1*p)>=0
    alpha = ZOOM(f,paramf,x,p,ff,g,a1,a0,c1,c2);
    return;
  end
  
  if i>maxit
    alpha = a1; % In this case, alpha satisfies only the first condition.
    fprintf('Fail to satisfy second Wolfe condition (outer loop)!\n');
    return;
  end
  
  a0 = a1;
  ff0 = ff1;
  a1 = rho*a1+(1-rho)*amax;
  i = i+1;
  
end


function alpha = ZOOM(f,paramf,x,p,ff,g,alo,ahi,c1,c2)
  
fflo = f((x+alo*p)',paramf{:});
maxit = 100;
j = 0;

while 1
  aj = (alo+ahi)/2;
  xj = x+aj*p;
  [ffj,gj] = f(xj',paramf{:});
  
  if (ffj>ff+c1*aj*g*p) | (ffj>=fflo)
    ahi = aj;
  else
    if abs(gj*p)<=-c2*g*p
      alpha = aj;
      return;
    end
    if gj*p*(ahi-alo)>=0
      ahi = alo;
    end
    alo = aj;
    fflo = ffj;
  end
  
  j = j + 1;
  if j==maxit
    alpha = aj; % In this case, alpha satisfies only the first condition.
    fprintf('Fail to satisfy second Wolfe condition (inner loop)!\n');
    return;
  end
end

