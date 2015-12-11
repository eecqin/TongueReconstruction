
addpath /home/cqin/matlab_files/src/rbf/

params.rbf.M = 500;
params.rbf.sigma = 55;
params.rbf.l = 1e-4;

% --- DATA 
%
% Load dataset
load feal0_adapt.mat
V2 = V; iTR2 = iTR; iTE2 = iTE;
load maaw0_adapt.mat
V = s1.V; iTR = s1.iTR; iTE = s1.iTE;

I = [2 9 19]';
II = [2*I-1 2*I]';                  % Indices of K landmarks
J = [1:24]'; JJ = [2*J-1 2*J]';     % Indices of P contour points
    
% Training and testing set
yTR = V(iTR,JJ(:)); xTR = yTR(:,II(:)); yTE = V(iTE,JJ(:)); xTE = yTE(:,II(:));
% Load predictive models
load(['f_K' num2str(length(I)) '_M' num2str(params.rbf.M) '_s' num2str(params.rbf.sigma) ...
          '_l' num2str(params.rbf.l) '.mat'],'f','fl','-mat');

% Test
hTEl = linf(xTE,fl);
hTE = rbf(xTE,f);

tonguermse(hTE,yTE),tonguermse(hTEl,yTE)
mean(tonguermse(hTE,yTE)),mean(tonguermse(hTEl,yTE))


