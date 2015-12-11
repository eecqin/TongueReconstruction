format compact


% --- PARAMETERS
params.adapt.iA0 = [0.35 -2.5; 0.5 0.35];
params.adapt.ib0 = [150; -50];
%params.adapt.iA0 = 0.8*[cos(pi/8) -sin(pi/8);sin(pi/8) cos(pi/8)];
%params.adapt.ib0 = [-3;-5]; 
params.adapt.tol = 1e-3;           % Set this carefully
params.rbf.M = 500;
params.rbf.sigma = 55;
params.rbf.l = 1e-4;
K = 3; P = 24;
% --- DATA 
load maaw0_adapt.mat
V = s1.V; iTR = s1.iTR; iTE = s1.iTE;
I = [2 9 19]'; J = [1:24]';
II = [2*I-1 2*I]'; JJ = [2*J-1 2*J]';
yTR = V(iTR,JJ(:)); xTR = yTR(:,II(:));
yTE = V(iTE,JJ(:)); xTE = yTE(:,II(:));

% iAb
blka = zeros(2*P,2*P);
ba = zeros(2*P,1);
iA0 = [0.35 -2.5; 0.5 0.35]; ib0 = [150; -50]; % original [5;-5]

iA1 = [1.1*iA0]; ib1 = [ib0+[10;-5]];
iA2 = [1.0*iA0]; ib2 = [ib0];
iA3 = [1.2*iA0]; ib3 = [ib0+[30;-10]];
for p=1:8
   blka(2*p-1:2*p,2*p-1:2*p) = iA1; % 
   ba(2*p-1:2*p) = ib1;
end
for p=9:16
   blka(2*p-1:2*p,2*p-1:2*p) = iA2; % 
   ba(2*p-1:2*p) = ib2;
end
for p=17:24
   blka(2*p-1:2*p,2*p-1:2*p) = iA3; % 
   ba(2*p-1:2*p) = ib3;
end

yTE2 = yTE*blka' + repmat(ba',size(yTE,1),1); xTE2 = yTE2(:,II(:));
figure(123);clf;
for i=1:size(yTE,1)
  set(0,'CurrentFigure',123); clf;
  plot(yTE(i,1:2:end-1),yTE(i,2:2:end),'bo-');hold on;
  plot(yTE2(i,1:2:end-1),yTE2(i,2:2:end),'ro-');
  set(gca,'DataAspectRatio',[1 1 1]); axis ij; 
  axis([-90 140 0 100]);
  legend('old speaker','new speaker','Location','NorthEast');
  drawnow;pause
end
% Adaptation and evaluation data
NN = 50;
Y_adpt2 = yTE2(1:NN,:); X_adpt2 = Y_adpt2(:,II(:)); 
Y_eval2 = yTE2(NN+1:end,:); X_eval2 = Y_eval2(:,II(:)); 

% --- MODEL
fl = linftrain(xTR,yTR);
load(['rbfcentres' num2str(params.rbf.M) '.mat'],'mu','-mat'); 
muI = mu(:,II(:)); % kmeans
f = rbftrain(xTR,yTR,muI,params.rbf.sigma,params.rbf.l);	 % Fit rbf  

return

% --- CORE
l = 0;
% fcn F for nonadaptive case by linf
tic; [Al,bl] = linfadapt(X_adpt2,Y_adpt2,fl,l); toc;
% fcn F for nonadaptive case rbf
tic; [A,b] = rbfadapt(X_adpt2,Y_adpt2,f,Al,bl,l,1e-6); toc; 

% initialization for adaptive cases
Abl = [repmat(Al(:),K,1); repmat(bl(:),K,1); repmat(Al(:),P,1); repmat(bl(:),P,1)]';
iAl = inv(Al); ibl = -inv(Al)*bl; Abl_2 = [repmat(Al(:),K,1); repmat(bl(:),K,1); repmat(iAl(:),P,1); repmat(ibl(:),P,1)]';
Ab = [repmat(A(:),K,1); repmat(b(:),K,1); repmat(A(:),P,1); repmat(b(:),P,1)]';
iA = inv(A); ib = -inv(A)*b; Ab_2 = [repmat(A(:),K,1); repmat(b(:),K,1); repmat(iA(:),P,1); repmat(ib(:),P,1)]';
% end of initialization for adaptive cases

% $$$ %%%%%%%%% fcn F for adaptive case by linf
% $$$ tic; [Abl2,El2] = linfadapt2(X_adpt2,Y_adpt2,fl,0); toc;
% $$$ %%%%%%%%% fcn F for adaptive case by rbf+BFGS with initialization from rbfadapt
% $$$ tic; [Ab2,E2] = rbfadapt2(X_adpt2,Y_adpt2,f,Ab',0,0.1); toc;

% fcn E for adaptive case by linf+BFGS with initialization from rbfadapt
tic; Abl2_2 = linfadapt22(X_adpt2,Y_adpt2,fl,Abl_2',l,1e-6); toc; % check
% fcn E for adaptive case by rbf+BFGS with initialization from rbfadapt
tic; Ab2_2 = rbfadapt22(X_adpt2,Y_adpt2,f,Ab_2',l,1e-6); toc;
% $$$ % fcn E for adaptive case by linf+AlterOpt with initialization from rbfadapt
% $$$ tic; Abl2_2a = rbfadapt22a(X_adpt2,Y_adpt2,fl,Ab_2',1e-6,500); toc;
% $$$ % fcn E for adaptive case by rbf+AlterOpt with initialization from rbfadapt
% $$$ tic; Ab2_2a = rbfadapt22a(X_adpt2,Y_adpt2,f,Ab_2',1e-6,500); toc;

% Evaluation set
[Yl,fgXl,gXl] = adapted_f(X_eval2,fl,Al,bl);              % fcn F for nonadaptive case by linf
[Y,fgX,gX] = adapted_f(X_eval2,f,A,b);                    % fcn F for nonadaptive case by rbf
[Yl2_2,fgXl2_2,gXl2_2] = adapted_f22(X_eval2,fl,Abl2_2); % fcn E for adaptive case by linf+BFGS
[Y2_2,fgX2_2,gX2_2] = adapted_f22(X_eval2,f,Ab2_2);      % fcn E for adaptive case by rbf+BFGS
% $$$ [Yl2_2a,fgXl2_2a,gXl2_2a] = adapted_f22(X_eval2,fl,Abl2_2a);  % fcn E for adaptive case by linf+AlterOpt
% $$$ [Y2_2a,fgX2_2a,gX2_2a] = adapted_f22(X_eval2,f,Ab2_2a);  % fcn E for adaptive case by rbf+AlterOpt


figure(123);clf;
for i=1:size(Y_eval2,1)
  set(0,'CurrentFigure',123);clf;
  subplot(321);
  plot(Y_eval2(i,1:2:end-1),Y_eval2(i,2:2:end),'bo-');hold on;
  plot(Yl(i,1:2:end-1),Yl(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
  axis([40 140 40 100]);
  title(['linfadapt on fcn F for nonadaptive case: ' num2str(norm(Yl-Y_eval2,'fro'))]);
  subplot(322);
  plot(Y_eval2(i,1:2:end-1),Y_eval2(i,2:2:end),'bo-');hold on;
  plot(Y(i,1:2:end-1),Y(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
  axis([40 140 40 100]);
  title(['rbfadapt (BFBS) on fcn F for nonadaptive case: ' num2str(norm(Y-Y_eval2,'fro'))]);
  subplot(323);
  plot(Y_eval2(i,1:2:end-1),Y_eval2(i,2:2:end),'bo-');hold on;
  plot(Yl2_2(i,1:2:end-1),Yl2_2(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
  axis([40 140 40 100]);
  title(['linfadapt (BFBS) on fcn E for adaptive case: ' num2str(norm(Yl2_2-Y_eval2,'fro'))]);
  subplot(324);
  plot(Y_eval2(i,1:2:end-1),Y_eval2(i,2:2:end),'bo-');hold on;
  plot(Y2_2(i,1:2:end-1),Y2_2(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
  axis([40 140 40 100]);
  title(['rbfadapt (BFGS) on fcn E for adaptive case: ' num2str(norm(Y2_2-Y_eval2,'fro'))]);
% $$$ subplot(325);
% $$$ plot(Y_eval2(i,1:2:end-1),Y_eval2(i,2:2:end),'bo-');hold on;
% $$$ plot(Yl2_2a(i,1:2:end-1),Yl2_2a(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
% $$$ axis([40 140 40 100]);
% $$$ title(['linfadapt (AlterOpt) on fcn E for adaptive case: ' num2str(norm(Yl2_2a-Y_eval2,'fro'))]);
% $$$ subplot(326);
% $$$ plot(Y_eval2(i,1:2:end-1),Y_eval2(i,2:2:end),'bo-');hold on;
% $$$ plot(Y2_2a(i,1:2:end-1),Y2_2a(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
% $$$ axis([40 140 40 100]);
% $$$ title(['rbfadapt (AlterOpt) on fcn E for adaptive case: ' num2str(norm(Y2_2a-Y_eval2,'fro'))]);
  drawnow; pause;
end
