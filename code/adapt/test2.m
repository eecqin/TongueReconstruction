format compact


% --- PARAMETERS
params.adapt.tol = 1e-6;           % Set this carefully
params.adapt.l = 1e+5;                % \lambda for regularization
params.rbf.M = 500;
params.rbf.sigma = 55;
params.rbf.l = 1e-2;
N = 10;

% --- DATA
load feal0_adapt.mat
V2 = V; iTR2 = iTR; iTE2 = iTE;
load maaw0_adapt.mat
V = s1.V; iTR = s1.iTR; iTE = s1.iTE;

% $$$ I = [3 12 21]';      % good
I = [2 9 19]';       % bad
% $$$ I = [2 8 14 20]';    % good
% $$$ I = [2 7 14 21]';    % bad
% $$$ I = [2 7 12 17 22]'; % good
% $$$ I = [2 6 12 19 23]'; % bad
J = [1:24]';
K = length(I); P = length(J);
II = [2*I-1 2*I]'; JJ = [2*J-1 2*J]';
yTR = V(iTR,JJ(:)); xTR = yTR(:,II(:)); yTE = V2(iTE2,JJ(:)); xTE = yTE(:,II(:));
NN = 500;  % NN contours for adaptation
Y_adpt0 = yTE(1:NN,:); X_adpt0 = Y_adpt0(:,II(:)); 
Y_eval0 = yTE(NN+1:end,:); X_eval0 = Y_eval0(:,II(:)); 
ii = randperm(NN);    % random permutation
Y_adpt2 = Y_adpt0(ii([1:N]),:); X_adpt2 = Y_adpt2(:,II(:));   % clean set to be transformed for adaptation
Y_eval2 = Y_eval0; X_eval2 = Y_eval2(:,II(:));

% sanity check
figure(1);clf;
subplot(121);plot(V(:,1:2:end-1),V(:,2:2:end),'b.');
title('tongue contours for maaw0');
set(gca,'DataAspectRatio',[1 1 1]);axis ij;axis([30 140 30 110]);
subplot(122);plot(V2(:,1:2:end-1),V2(:,2:2:end),'r.');
title('tongue contours for feal0')
set(gca,'DataAspectRatio',[1 1 1]);axis ij;axis([30 140 30 110]);
figure(2);clf;
plot(Y_eval2(:,1:2:end-1),Y_eval2(:,2:2:end),'b.');hold on;
plot(Y_adpt2(:,1:2:end-1),Y_adpt2(:,2:2:end),'r.');
title('adaptation set again evaluation set both from feal0');
set(gca,'DataAspectRatio',[1 1 1]);axis ij;
figure(3);clf;
plot(xTR(:,1:2:end-1),xTR(:,2:2:end),'b.');hold on;
plot(X_eval2(:,1:2:end-1),X_eval2(:,2:2:end),'r.');
title('Scatterplot of landmarks from old and new speaker');
set(gca,'DataAspectRatio',[1 1 1]);axis ij;
% end of sanity check


% --- MODEL
fl = linftrain(xTR,yTR);
load(['rbfcentres' num2str(params.rbf.M) '.mat'],'mu','-mat'); 
muI = mu(:,II(:)); % kmeans
f = rbftrain(xTR,yTR,muI,params.rbf.sigma,params.rbf.l);	 % Fit rbf  


% --- CORE
% fcn F for nonadaptive case by linf
tic; [Al,bl] = linfadapt(X_adpt2,Y_adpt2,fl,params.adapt.l); toc;
% $$$ % fcn F for nonadaptive case rbf
tic; [A,b] = rbfadapt(X_adpt2,Y_adpt2,f,Al,bl,params.adapt.l,params.adapt.tol); toc; 

% initialization for adaptive cases
Abl = [repmat(Al(:),K,1); repmat(bl(:),K,1); repmat(Al(:),P,1); repmat(bl(:),P,1)]';
iAl = inv(Al); ibl = -inv(Al)*bl; Abl_2 = [repmat(Al(:),K,1); repmat(bl(:),K,1); repmat(iAl(:),P,1); repmat(ibl(:),P,1)]';
Ab = [repmat(A(:),K,1); repmat(b(:),K,1); repmat(A(:),P,1); repmat(b(:),P,1)]';
iA = inv(A); ib = -inv(A)*b; Ab_2 = [repmat(A(:),K,1); repmat(b(:),K,1); repmat(iA(:),P,1); repmat(ib(:),P,1)]';
% end of initialization for adaptive cases

% fcn E for adaptive case by rbf+BFGS with initialization from rbfadapt
tic; Ab2_2 = rbfadapt22(X_adpt2,Y_adpt2,f,Ab_2',params.adapt.l,params.adapt.tol); toc;


[Yl,fgXl,gXl] = adapted_f(X_eval2,fl,Al,bl);              % fcn F for nonadaptive case by linf
[Y,fgX,gX] = adapted_f(X_eval2,f,A,b);                    % fcn F for nonadaptive case by rbf
[Y2_2,fgX2_2,gX2_2] = adapted_f22(X_eval2,f,Ab2_2);      % fcn E for adaptive case by rbf+BFGS
mean(tonguermse(Yl,Y_eval2)),mean(tonguermse(Y,Y_eval2)),mean(tonguermse(Y2_2,Y_eval2))

i = 100;figure(123);clf;
subplot(221);
plot(Y_eval2(i,1:2:end-1),Y_eval2(i,2:2:end),'bo-');hold on;
plot(Yl(i,1:2:end-1),Yl(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
axis([40 140 40 100]);
title(['linfadapt on fcn F for nonadaptive case: ' num2str(norm(Yl-Y_eval2,'fro'))]);
subplot(222);
plot(Y_eval2(i,1:2:end-1),Y_eval2(i,2:2:end),'bo-');hold on;
plot(Y(i,1:2:end-1),Y(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
axis([40 140 40 100]);
title(['rbfadapt (BFBS) on fcn F for nonadaptive case: ' num2str(norm(Y-Y_eval2,'fro'))]);
subplot(224);
plot(Y_eval2(i,1:2:end-1),Y_eval2(i,2:2:end),'bo-');hold on;
plot(Y2_2(i,1:2:end-1),Y2_2(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
axis([40 140 40 100]);
title(['rbfadapt (BFBS) on fcn E for adaptive case: ' num2str(norm(Y2_2-Y_eval2,'fro'))]);

return

% fcn E for adaptive case by linf+BFGS with initialization from rbfadapt
tic; Abl2_2 = linfadapt22(X_adpt2,Y_adpt2,fl,Abl_2',params.adapt.l,params.adapt.tol); toc;
% fcn E for adaptive case by rbf+BFGS with initialization from rbfadapt
tic; Ab2_2 = rbfadapt22(X_adpt2,Y_adpt2,f,Ab_2',params.adapt.l,params.adapt.tol); toc;

% Evaluation set
[Yl,fgXl,gXl] = adapted_f(X_eval2,fl,Al,bl);              % fcn F for nonadaptive case by linf
[Y,fgX,gX] = adapted_f(X_eval2,f,A,b);                    % fcn F for nonadaptive case by rbf
[Yl2_2,fgXl2_2,gXl2_2] = adapted_f22(X_eval2,fl,Abl2_2); % fcn E for adaptive case by linf+BFGS
[Y2_2,fgX2_2,gX2_2] = adapted_f22(X_eval2,f,Ab2_2);      % fcn E for adaptive case by rbf+BFGS



figure(123);clf;
for i=1:size(Y_eval2,1)
  set(0,'CurrentFigure',123);clf;
  subplot(221);
  plot(Y_eval2(i,1:2:end-1),Y_eval2(i,2:2:end),'bo-');hold on;
  plot(Yl(i,1:2:end-1),Yl(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
  axis([40 140 40 100]);
  title(['linfadapt on fcn F for nonadaptive case: ' num2str(norm(Yl-Y_eval2,'fro'))]);
  subplot(222);
  plot(Y_eval2(i,1:2:end-1),Y_eval2(i,2:2:end),'bo-');hold on;
  plot(Y(i,1:2:end-1),Y(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
  axis([40 140 40 100]);
  title(['rbfadapt (BFBS) on fcn F for nonadaptive case: ' num2str(norm(Y-Y_eval2,'fro'))]);
% $$$   subplot(223);
% $$$   plot(Y_eval2(i,1:2:end-1),Y_eval2(i,2:2:end),'bo-');hold on;
% $$$   plot(Yl2_2(i,1:2:end-1),Yl2_2(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
% $$$   axis([40 140 40 100]);
% $$$   title(['linfadapt (BFBS) on fcn E for adaptive case: ' num2str(norm(Yl2_2-Y_eval2,'fro'))]);
% $$$   subplot(224);
% $$$   plot(Y_eval2(i,1:2:end-1),Y_eval2(i,2:2:end),'bo-');hold on;
% $$$   plot(Y2_2(i,1:2:end-1),Y2_2(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
% $$$   axis([40 140 40 100]);
% $$$   title(['rbfadapt (BFGS) on fcn E for adaptive case: ' num2str(norm(Y2_2-Y_eval2,'fro'))]);
  drawnow; pause;
end
