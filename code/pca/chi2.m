% chi2(y,n) Area under a chi-squared distribution
%
% Computes the area under a chi-squared distribution with n degrees of
% freedom for 0<x<y, i.e.
%     chi2(y,n) = int^y_0{f(x;n) dx}
% where
%     f(x;n) = exp((n/2-1)*log(x)-x/2-n*log(2)/2-lngamma(n/2))
% is the density of the chi-squared distribution with n degrees of freedom.

function c = chi2(y,n)

c = gammainc(y/2,n/2);
