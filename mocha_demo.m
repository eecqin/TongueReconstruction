addpath(genpath('code'));

% Load tongue contour set
disp('load tongue contours');
load dat/tongue.mat;

% Fit rbf model by less optimal but fast training
disp('fit kmeans to all tongue contour points');
nTR = 5000; yTR = Y(1:nTR,:);
[mu,labels,errfunc,code] = kmeans(yTR,500,'rndmeans',100,1e-6);
disp('fit the model for the landmarks of interest');
I = [4 9 14]'; II = [2*I-1 2*I]';
muI = mu(:,II(:)); xTR = yTR(:,II(:));
f = rbftrain(xTR,yTR,muI,55,1e-2); f0 = f; f.W = f.W(II(:),:);

% Adapt the rbf model
disp('load adaptation data for a new speaker');
load('dat/fsew0.mat','ema');
numAdapt=5000; X_adpt2 = ema(1:numAdapt,[7:12]);
lambda=1e+4; tol = 1e-6; A0 = eye(2); b0 = zeros(2,1);    % zero initial since no local optima
disp('adapt the rbf model');
tic; [A,b,E] = rbfadapt(X_adpt2,X_adpt2,f,A0,b0,lambda,tol); toc;
A,b,cond(A)

% Apply the adapted model to reconstruct tongue contours for the new speaker
X_eval2 = ema(numAdapt+1:numAdapt+50000,[7:12]); Y_eval2 = X_eval2;
[fX_eval2,fX_eval1,X_eval1] = adapted_f(X_eval2, f0, A, b);

% Load the hard constraints, e.g. pal and pha
pal=load(['dat/fsew0_pal.dat'],'-ascii');

% Plot tongue reconstructions for randomly sampled frames
figure(1); plotrec(fX_eval2,Y_eval2,[],[],pal,'MOCHA');




