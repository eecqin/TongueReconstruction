format compact


% --- PARAMETERS
params.adapt.tol = 1e-6;           % Set this carefully
params.rbf.M = 500;
params.rbf.sigma = 55;
params.rbf.l = 1e-4;
K = 3; P = 24;
% --- DATA 
load maaw0_adapt.mat
V = s1.V; iTR = s1.iTR; iTE = s1.iTE;

I = [2 9 19]'; J = [1:24]';
II = [2*I-1 2*I]'; JJ = [2*J-1 2*J]';
yTR = V(iTR,JJ(:)); xTR = yTR(:,II(:)); yTE = V(iTE,JJ(:)); xTE = yTE(:,II(:));
Y_adpt2 = yTE(1:50,:); X_adpt2 = Y_adpt2(:,II(:)); 
Y_eval2 = yTE(51:end,:); X_eval2 = Y_eval2(:,II(:)); 

% transform testing set
% $$$ % centre and scale dataset
% $$$ [UUU,VVV,mmm] = pca(yTE,size(yTE,2));
% $$$ figure(123);clf
% $$$ for i=1:48
% $$$   plot(UUU(1:2:end-1,i),UUU(2:2:end,i),'bo-');hold on;
% $$$   set(gca,'DataAspectRatio',[1 1 1]);pause
% $$$ end
[mu,labels,errfunc,code] = kmeans(yTE,20,[],100,1e-5);
[mu,labels,errfunc,code] = kmeans(yTE,50);
figure(2);clf;
% $$$ plot(mu(:,1:2:end-1),mu(:,2:2:end),'m.');
yTE = yTE(randperm(size(yTE,1)),:);
for i=10:550
  %w = randn(1,50);
  %w = 1 + 1*randn(1,50); %w = w./sum(w);
  %mu2 = w*mu;
  %plot(mu2(1:2:end-1),mu2(2:2:end),'bo-','LineWidth',3);
  plot(mean(yTE(i-2:i+2,1:2:end-1),1),mean(yTE(i-2:i+2,2:2:end),1),'bo-','LineWidth',1);
  set(gca,'DAtaAspectRatio',[1 1 1]); axis ij
  drawnow;pause
end

figure
plot(yTR(:,1:2:end-1),yTR(:,2:2:end),'b.');hold on;
plot(yTE(:,1:2:end-1),yTE(:,2:2:end),'r.');
set(gca,'DAtaAspectRatio',[1 1 1]); axis ij

% --- MODEL
fl = linftrain(xTR,yTR);
load(['rbfcentres' num2str(params.rbf.M) '.mat'],'mu','-mat'); 
muI = mu(:,II(:)); % kmeans
f = rbftrain(xTR,yTR,muI,params.rbf.sigma,params.rbf.l);	 % Fit rbf  

return

% --- CORE
NN = 50; l = 0;
% fcn F for nonadaptive case by linf
tic; [Al,bl] = linfadapt(X_adpt2(1:NN,:),Y_adpt2(1:NN,:),fl,l); toc;
% fcn F for nonadaptive case rbf
tic; [A,b] = rbfadapt(X_adpt2(1:NN,:),Y_adpt2(1:NN,:),f,Al,bl,l,1e-6); toc; 

% initialization for adaptive cases
Abl = [repmat(Al(:),K,1); repmat(bl(:),K,1); repmat(Al(:),P,1); repmat(bl(:),P,1)]';
iAl = inv(Al); ibl = -inv(Al)*bl; Abl_2 = [repmat(Al(:),K,1); repmat(bl(:),K,1); repmat(iAl(:),P,1); repmat(ibl(:),P,1)]';
Ab = [repmat(A(:),K,1); repmat(b(:),K,1); repmat(A(:),P,1); repmat(b(:),P,1)]';
iA = inv(A); ib = -inv(A)*b; Ab_2 = [repmat(A(:),K,1); repmat(b(:),K,1); repmat(iA(:),P,1); repmat(ib(:),P,1)]';
% end of initialization for adaptive cases

% $$$ %%%%%%%%% fcn F for adaptive case by linf
% $$$ tic; [Abl2,El2] = linfadapt2(X_adpt2(1:NN,:),Y_adpt2(1:NN,:),fl,0); toc;
% $$$ %%%%%%%%% fcn F for adaptive case by rbf+BFGS with initialization from rbfadapt
% $$$ tic; [Ab2,E2] = rbfadapt2(X_adpt2(1:NN,:),Y_adpt2(1:NN,:),f,Ab',0,0.1); toc;

% fcn E for adaptive case by linf+BFGS with initialization from rbfadapt
tic; Abl2_2 = linfadapt22(X_adpt2(1:NN,:),Y_adpt2(1:NN,:),fl,Abl_2',l,1e-6); toc;
% fcn E for adaptive case by rbf+BFGS with initialization from rbfadapt
tic; Ab2_2 = rbfadapt22(X_adpt2(1:NN,:),Y_adpt2(1:NN,:),f,Ab_2',l,1e-6); toc;
% $$$ % fcn E for adaptive case by linf+AlterOpt with initialization from rbfadapt
% $$$ tic; Abl2_2a = rbfadapt22a(X_adpt2(1:NN,:),Y_adpt2(1:NN,:),fl,Ab_2',1e-6,500); toc;
% $$$ % fcn E for adaptive case by rbf+AlterOpt with initialization from rbfadapt
% $$$ tic; Ab2_2a = rbfadapt22a(X_adpt2(1:NN,:),Y_adpt2(1:NN,:),f,Ab_2',1e-6,500); toc;

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
  title(['linfadapt on fcn F for nonadaptive case: ' num2str(norm(Yl-Y_eval2))]);
  subplot(322);
  plot(Y_eval2(i,1:2:end-1),Y_eval2(i,2:2:end),'bo-');hold on;
  plot(Y(i,1:2:end-1),Y(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
  axis([40 140 40 100]);
  title(['rbfadapt (BFBS) on fcn F for nonadaptive case: ' num2str(norm(Y-Y_eval2))]);
  subplot(323);
  plot(Y_eval2(i,1:2:end-1),Y_eval2(i,2:2:end),'bo-');hold on;
  plot(Yl2_2(i,1:2:end-1),Yl2_2(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
  axis([40 140 40 100]);
  title(['linfadapt (BFBS) on fcn E for adaptive case: ' num2str(norm(Yl2_2-Y_eval2))]);
  subplot(324);
  plot(Y_eval2(i,1:2:end-1),Y_eval2(i,2:2:end),'bo-');hold on;
  plot(Y2_2(i,1:2:end-1),Y2_2(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
  axis([40 140 40 100]);
  title(['rbfadapt (BFGS) on fcn E for adaptive case: ' num2str(norm(Y2_2-Y_eval2))]);
% $$$ subplot(325);
% $$$ plot(Y_eval2(i,1:2:end-1),Y_eval2(i,2:2:end),'bo-');hold on;
% $$$ plot(Yl2_2a(i,1:2:end-1),Yl2_2a(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
% $$$ axis([40 140 40 100]);
% $$$ title(['linfadapt (AlterOpt) on fcn E for adaptive case: ' num2str(norm(Yl2_2a-Y_eval2))]);
% $$$ subplot(326);
% $$$ plot(Y_eval2(i,1:2:end-1),Y_eval2(i,2:2:end),'bo-');hold on;
% $$$ plot(Y2_2a(i,1:2:end-1),Y2_2a(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
% $$$ axis([40 140 40 100]);
% $$$ title(['rbfadapt (AlterOpt) on fcn E for adaptive case: ' num2str(norm(Y2_2a-Y_eval2))]);
  drawnow; pause;
end

return

figure(125);clf;
subplot(121);plot(Abl_2,'b-');hold on;plot(Abl2_2,'r-');plot(Abl2_2a,'m-');
title('Recovered transformations by linfadapt'); xlabel('Dimension'); ylabel('Magnitude');
legend('linfadapt (BFBS )on fcn F for nonadaptive case',...
       'linfadapt (BFGS) on fcn E for adaptive case',...
       'linfadapt (AlterOpt) on fcn E for adaptive case','Location','NorthEast');
subplot(122);plot(Ab_2,'b-');hold on;plot(Ab2_2,'r-');plot(Ab2_2a,'m-');
title('Recovered transformations by rbfadapt'); xlabel('Dimension'); ylabel('Magnitude');
legend('rbfadapt (BFBS )on fcn F for nonadaptive case',...
       'rbfadapt (BFGS) on fcn E for adaptive case',...
       'rbfadapt (AlterOpt) on fcn E for adaptive case','Location','NorthEast');


% check on condition number
% rbf case
cond_Ab = 0; Ax_Ab = Ab(1:4*K); Cy_Ab = Ab(6*K+1:6*K+4*P);
cond_Ab2_2 = 0; Ax_Ab2_2 = Ab2_2(1:4*K); Cy_Ab2_2 = Ab2_2(6*K+1:6*K+4*P);
cond_Ab2_2a = 0; Ax_Ab2_2a = Ab2_2a(1:4*K); Cy_Ab2_2a = Ab2_2a(6*K+1:6*K+4*P);
for k=1:K
  cond_Ab = cond_Ab + cond(reshape(Ax_Ab(4*(k-1)+1:4*k),2,2));
  cond_Ab2_2 = cond_Ab2_2 + cond(reshape(Ax_Ab2_2(4*(k-1)+1:4*k),2,2));
  cond_Ab2_2a = cond_Ab2_2a + cond(reshape(Ax_Ab2_2a(4*(k-1)+1:4*k),2,2));
end
for p=1:P
  cond_Ab = cond_Ab + cond(reshape(Cy_Ab(4*(p-1)+1:4*p),2,2));
  cond_Ab2_2 = cond_Ab2_2 + cond(reshape(Cy_Ab2_2(4*(p-1)+1:4*p),2,2));
  cond_Ab2_2a = cond_Ab2_2a + cond(reshape(Cy_Ab2_2a(4*(p-1)+1:4*p),2,2));
end
cond_Ab/(K+P),cond_Ab2_2/(K+P),cond_Ab2_2a/(K+P)
% linf case
cond_Abl = 0; Ax_Abl = Abl(1:4*K); Cy_Abl = Abl(6*K+1:6*K+4*P);
cond_Abl2_2 = 0; Ax_Abl2_2 = Abl2_2(1:4*K); Cy_Abl2_2 = Abl2_2(6*K+1:6*K+4*P);
cond_Abl2_2a = 0; Ax_Abl2_2a = Abl2_2a(1:4*K); Cy_Abl2_2a = Abl2_2a(6*K+1:6*K+4*P);
for k=1:K
  cond_Abl = cond_Abl + cond(reshape(Ax_Abl(4*(k-1)+1:4*k),2,2));
  cond_Abl2_2 = cond_Abl2_2 + cond(reshape(Ax_Abl2_2(4*(k-1)+1:4*k),2,2));
  cond_Abl2_2a = cond_Abl2_2a + cond(reshape(Ax_Abl2_2a(4*(k-1)+1:4*k),2,2));
end
for p=1:P
  cond_Abl = cond_Abl + cond(reshape(Cy_Abl(4*(p-1)+1:4*p),2,2));
  cond_Abl2_2 = cond_Abl2_2 + cond(reshape(Cy_Abl2_2(4*(p-1)+1:4*p),2,2));
  cond_Abl2_2a = cond_Abl2_2a + cond(reshape(Cy_Abl2_2a(4*(p-1)+1:4*p),2,2));
end
cond_Abl/(K+P),cond_Abl2_2/(K+P),cond_Abl2_2a/(K+P)
