% [y,xx]=chi2inv(x,n) Inverse of chi2(y,n)
%
% By the bisection method, finds y = the inverse of chi2(x,n).
% Only works for x in [0.7,0.95]. xx = chi2(y,n) should be approximately x.

function [y,xx] = chi2inv(x,n)

% [a,b] is heuristically computed so that the sought value is in [a,b]
% for x in [0.7,0.95]. For x>0.95 b must be larger; for x<0.7 a must
% be smaller.

a=n;
b=a+2.1*sqrt(2*n);
c=0;
while abs(chi2(c,n)-x)>=0.00005
  c=(a+b)/2;
  if chi2(c,n)>x
    b=c;
  else
    a=c;
  end
end

y=c;

if nargout == 2
  xx=chi2(y,n);
end

