% [y,F] = XRMBtrans(y_old) dynamic model in articulatory data

function [y,F] = XRMBtrans(y_old)

% Mixture of linear models
% $$$ load flinf.mat
% $$$ y = linf(y_old,f); 
% $$$ F = f.W;
% $$$ % $$$ % Alternatively
% $$$ % $$$ y = y_old;
% $$$ % $$$ if (nargout>1) F = eye(16); end
  
% RBF
load flinf.mat
[y,F] = linf(y_old,f);
