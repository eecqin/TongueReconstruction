addpath(genpath('code'));

% Load tongue contour set
disp('load tongue contours');
load dat/tongue.mat;

% Fit rbf model by less optimal but fast training
disp('fit kmeans to all tongue contour points');
nTR = 5000; yTR = Y(1:nTR,:);
[mu,labels,errfunc,code] = kmeans(yTR,500,'rndmeans',100,1e-5);
disp('fit the model for the landmarks of interest');
I = [3 7 11 15]'; II = [2*I-1 2*I]';
muI = mu(:,II(:)); xTR = yTR(:,II(:));
f = rbftrain(xTR,yTR,muI,55,1e-1); f0 = f; f.W = f.W(II(:),:);

% Adapt the rbf model
disp('load adaptation data for a new speaker');
load('dat/jw11.mat','xrmb');
numAdapt=5000; X_adpt2 = xrmb(1:numAdapt,[5:12]);
lambda=1e+5; tol = 1e-6; A0 = eye(2); b0 = zeros(2,1);    % zero initial since no local optima
disp('adapt the rbf model');
tic; [A,b,E] = rbfadapt(X_adpt2,X_adpt2,f,A0,b0,lambda,tol); toc;
A,b,cond(A)

% Apply the adapted model to reconstruct tongue contours for the new speaker
X_eval2 = xrmb(numAdapt+1:numAdapt+20000,[5:12]); Y_eval2 = X_eval2;
[fX_eval2,fX_eval1,X_eval1] = adapted_f(X_eval2, f0, A, b);

% Load the hard constraints, e.g. pal and pha
pal=load(['dat/jw11_pal.dat'],'-ascii'); pal=pal/1000;
pha=load(['dat/jw11_pha.dat'],'-ascii'); pha=pha/1000;

% Plot tongue reconstructions for randomly sampled frames
figure(1); plotrec(fX_eval2,Y_eval2,[],pha,pal,'XRMB');


