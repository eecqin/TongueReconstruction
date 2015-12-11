function [xy,curvexy,kk] = my_cspline(dd,I,nP,min_x,max_x)

kk = dd(:,I);

curvexy = zeros(nP,2); 
curvexy(:,1) = linspace(min_x,max_x,nP);
%curvexy = dd';  % ground truth for debugging
%curvexy(:,2) = interp1(kk(1,:),kk(2,:),curvexy(:,1)','pchip','extrap')';
%curvexy(:,2) = interp1(kk(1,:),kk(2,:),curvexy(:,1)','spline','extrap')';
curvexy(:,2) = interp1(kk(1,:),kk(2,:),curvexy(:,1)','spline')';

% distance2curve (very slow)
%xy = distance2curve(curvexy,dd',setdiff([1:24],I),'pchip');
%xy(I',:) = dd(:,I)';

% ALTERNATIVE DISCRETIZED VERSION IMPLEMENTATION
sqd = sqdist(curvexy,dd');
[tmp,II] = min(sqd,[],1);
xy = curvexy(II,:);
