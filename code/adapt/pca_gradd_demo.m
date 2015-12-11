% Script to demo adaptation using PCA as front-end to rescale adaptation set 
% (and evaluation set) and then apply gradient descent algorithm to minimize 
% the cost function. And it works!

clear all

format compact
% ----------------------- Data ------------------------------------------------- %
% Load dataset
load maaw0_adapt.mat
V = s1.V; iTR = s1.iTR; iTE = s1.iTE;

% $$$ % Plot original data
% $$$ figure(1);clf;
% $$$ plot(V(iTR,[1:2:47]),V(iTR,[2:2:48]),'b.',V(iTE,[1:2:47]),V(iTE,[2:2:48]),'r.');
% $$$ axis([20 140 20 120]); axis('ij'); set(gca,'DataAspectRatio',[1 1 1]);

% ----------------------------- Train and test RBF ----------------------------- %
% Indices of contour points
I = [2 9 19]'; J = [1:24]';
II = [2*I-1 2*I]'; JJ = [2*J-1 2*J]';
yTR = V(iTR,JJ(:)); xTR = yTR(:,II(:)); yTE = V(iTE,JJ(:)); xTE = yTE(:,II(:));

% Centre and scale the training set
XXX = yTR(:,[1:2:end-1]); YYY = yTR(:,[2:2:end]);
[UUU,VVV,mmm] = pca([XXX(:) YYY(:)],2);
yTR = (yTR-repmat(mmm,size(yTR,1),P))./max(sqrt(VVV'));
xTR = yTR(:,II(:));
yTEn = (yTE-repmat(mmm,size(yTE,1),P))./max(sqrt(VVV'));
xTEn = yTEn(:,II(:));

% Train RBF 
M = 500; sigma = 55; l = 1e-4;
% linf
fl = linftrain(xTR,yTR);
hlTEn = rbf(xTEn,fl);
hTEn = rbf(xTEn,f);
hTE = hTEn.*repmat(sqrt(VVV'),size(hTEn,1),P)+repmat(mmm,size(hTEn,1),P);
norm(hTE-yTE,'fro')

% rbf
%load(['rbfcentres' num2str(M) '.mat'],'mu','-mat'); 
%muI = mu(:,II(:)); % kmeans
f = rbftrain(xTR,yTR,M,sigma,l);	 % Fit rbf  
hTEn = rbf(xTEn,f);
hTE = hTEn.*repmat(sqrt(VVV'),size(hTEn,1),P)+repmat(mmm,size(hTEn,1),P);
norm(hTE-yTE,'fro')

return

figure(2);
for i=1:10:length(iTE)
  set(0,'CurrentFigure',2);clf;  
  plot(yTE(i,[1:2:47]),yTE(i,[2:2:48]),'co-','LineWidth',2,'MarkerFaceColor','c');hold on;
  plot(hTE(i,[1:2:47]),hTE(i,[2:2:48]),'ro-','LineWidth',2,'MarkerFaceColor','r');
  plot(xTE(i,[1:2:5]),xTE(i,[2:2:6]),'bo','LineWidth',2,'MarkerFaceColor','b');hold off;
  set(gca,'DataAspectRatio',[1 1 1]); axis([20 140 20 120]); axis ij; drawnow; pause
end

% ----------------------------------------------------------------------------- % 

% ----------------------- Adaptation ------------------------------------------ % 
tol1  = 1e-4;                         % precision on objective function
tol2  = 1e-4;                         % precision on parameters
eta   = 1e-9;                         % stepsize
maxit = 2000;                         % max num iterations
N = 10;                                % num. adaptation data

% Adaptation set (original)
Y_adpt = yTE([1:N]',:); X_adpt = Y_adpt(:,II(:));
% Independent evaluation set (original)
Y_eval = yTE([201:end]',:); X_eval = Y_eval(:,II(:));

MU = f.C; WW = f.W(:,1:end-1)'; vv = f.W(:,end)'; ss = f.s;
K = size(X_adpt,2)/2; P = size(Y_adpt,2)/2; M = size(WW,1); 

% 2D-wise transformation
iA0 = [2 1.5; -1.5 1]; ib0 = [50; 50];  % worked examples
% $$$ iA0 = [10 -1.5; 1.5 6]; ib0 = [50; 50]; 
% $$$ iA0 = [15 -10; -1.5 5]; ib0 = [100; -50]; 
A0 = inv(iA0); b0 = -inv(iA0)*ib0; 
Y_adpt2 = (kron(eye(P),iA0)*Y_adpt')'+repmat(repmat(ib0',1,P),size(Y_adpt,1),1); X_adpt2 = Y_adpt2(:,II(:));
Y_eval2 = (kron(eye(P),iA0)*Y_eval')'+repmat(repmat(ib0',1,P),size(Y_eval,1),1); X_eval2 = Y_eval2(:,II(:));

% Centre and scale dataset
XXX = Y_adpt2(:,[1:2:end-1]); YYY = Y_adpt2(:,[2:2:end]);
[UUU,VVV,mmm] = pca([XXX(:) YYY(:)],2);
%Y_adpt3 = (Y_adpt2-repmat(mmm,size(Y_adpt,1),P))./repmat(sqrt(VVV'),size(Y_adpt,1),P);
Y_adpt3 = (Y_adpt2-repmat(mmm,size(Y_adpt,1),P))./max(sqrt(VVV'));
X_adpt3 = Y_adpt3(:,II(:));
%Y_eval3 = (Y_eval2-repmat(mmm,size(Y_eval,1),P))./repmat(sqrt(VVV'),size(Y_eval,1),P);
Y_eval3 = (Y_eval2-repmat(mmm,size(Y_eval,1),P))./max(sqrt(VVV'));
X_eval3 = Y_eval3(:,II(:));

Y_adpt3 = Y_adpt2; X_adpt3 = X_adpt2;

% SANITY CHECK
figure(3);clf;
plot(Y_adpt(1,[1:2:end-1]),Y_adpt(1,[2:2:end]),'b.',...       % original
     Y_adpt2(1,[1:2:end-1]),Y_adpt2(1,[2:2:end]),'r.',...     % transformed db
     Y_adpt3(1,[1:2:end-1]),Y_adpt3(1,[2:2:end]),'m.',...     % transformed and scaled db
     Y_adpt(:,[1:2:end-1]),Y_adpt(:,[2:2:end]),'b.',...       
     Y_adpt2(:,[1:2:end-1]),Y_adpt2(:,[2:2:end]),'r.',...     
     Y_adpt3(:,[1:2:end-1]),Y_adpt3(:,[2:2:end]),'m.');
set(gca,'Box','on','DataAspectRatio',[1 1 1]);
legend('original','transformed','rescaled','Location','NorthEast'); axis('ij'); drawnow;

%%%%%%%%%%%%%%%%%%%%%  Adaptation using Gradient Descent %%%%%%%%%%%%%%%%%%
% Initial solution (doesn't matter)
A = eye(2); b = zeros(2,1);         % Neutral initial guess
% $$$ A = 100*randn(2,2); b = 2000*randn(2,1);   % Random initial 

old_A = 100*eye(2); old_b = [100 100]';
old_F = Fval(X_adpt3,Y_adpt3,WW,vv,MU,ss,old_A,old_b)/N;
F = Fval(X_adpt3,Y_adpt3,WW,vv,MU,ss,A,b)/N;
it = 1; 
while abs(F-old_F)>tol1 & norm([A b]-[old_A old_b])>tol2 & it<=maxit
  old_F = F; old_A = A; old_b = b;
  
  % Gradient descent iterations
  [gA,gb] = Fgrad(X_adpt3,Y_adpt3,WW,vv,MU,ss,A,b);
  
  A = A - eta*reshape(gA,2,2);
  b = b - eta*gb';
  F = Fval(X_adpt3,Y_adpt3,WW,vv,MU,ss,A,b)/N;
  
  % Echo the results
  fprintf(1,'Cycle %d: F=%1.6f, normD=%1.6f, A=[%1.3f %1.3f;%1.3f %1.3f], b=[%1.3f %1.3f]''\n',...
          it,F,norm([A b]-[old_A old_b]),A(1,1),A(1,2),A(2,1),A(2,2),b(1),b(2));
  
  % SANITY CHECK
% $$$   if mod(it,100)==0
% $$$     % Recovered dataset
% $$$     YY_adpt = (kron(eye(P),A)*Y_adpt3' + repmat(kron(ones(P,1),b),1,N))';
% $$$     XX_adpt = YY_adpt(:,II(:));
% $$$     % Plot
% $$$     figure(4);clf;
% $$$     plot(Y_adpt(1,[1:2:end-1]),Y_adpt(1,[2:2:end]),'b.',...  % original
% $$$          YY_adpt(1,[1:2:end-1]),YY_adpt(1,[2:2:end]),'m.',...          % recovered
% $$$          Y_adpt(:,[1:2:end-1]),Y_adpt(:,[2:2:end]),'b.',...
% $$$          YY_adpt(:,[1:2:end-1]),YY_adpt(:,[2:2:end]),'m.');
% $$$     set(gca,'Box','on','DataAspectRatio',[1 1 1]);
% $$$     legend('original','recovered','Location','NorthEast')
% $$$     axis('ij'); drawnow; pause;
% $$$   end
  it = it + 1;
end
% ------------------------------------------------------------------------------ %

% ---------------- Adaptation using BFGS ----------------------------- %
A = eye(2); b = zeros(2,1);         % Neutral initial guess
tol = 0.01;
[A,b] = rbfadapt(X_adpt2,Y_adpt2,f,A,b,0,tol);
% -------------------------------------------------------------------- %

% To use the adapted model
XX_adpt = (kron(eye(K),A)*X_adpt3' + repmat(kron(ones(K,1),b),1,size(X_adpt,1)))'; YY_adpt = rbf(XX_adpt,f);
XX_eval = (kron(eye(K),A)*X_eval3' + repmat(kron(ones(K,1),b),1,size(X_eval,1)))'; YY_eval = rbf(XX_eval,f);

% Transform back
YY_adpt3 = (kron(eye(P),inv(A))*YY_adpt' + repmat(kron(ones(P,1),-inv(A)*b),1,size(Y_adpt,1)))';
YY_eval3 = (kron(eye(P),inv(A))*YY_eval' + repmat(kron(ones(P,1),-inv(A)*b),1,size(Y_eval,1)))';

% Scale back
fX_adpt2 = YY_adpt3.*repmat(sqrt(VVV'),size(Y_adpt,1),P)+repmat(mmm,size(Y_adpt,1),P);
fX_eval2 = YY_eval3.*repmat(sqrt(VVV'),size(Y_eval,1),P)+repmat(mmm,size(Y_eval,1),P);

% compute RMSE 
E_adpt2 = tonguermse(fX_adpt2,Y_adpt2); E_eval2 = tonguermse(fX_eval2,Y_eval2);
fprintf(1,'e_adpt2=%1.3fmm e_eval2=%1.3fmm\n',mean(E_adpt2),mean(E_eval2)) 


% Plot results
figure(5);clf;
for i=1:N
  subplot(1,2,1);
  plot(Y_adpt2(i,1:2:end),Y_adpt2(i,2:2:end),'bo-',fX_adpt2(i,1:2:end),fX_adpt2(i,2:2:end),'ro-');   
  title(['rbf: i=' num2str(i) '; e2=' num2str(mean(tonguermse(fX_adpt2(i,:),Y_adpt2(i,:))),'%1.3fmm')]); 
  axis([min(min(Y_adpt2(:,1:2:end))) max(max(Y_adpt2(:,1:2:end))) min(min(Y_adpt2(:,2:2:end))) max(max(Y_adpt2(:,2:2:end)))]); 
  set(gca,'DataAspectRatio',[1 1 1]); axis ij; drawnow; 
  pause(0.5);
end
for i=1:50:size(Y_eval2,1)
  subplot(1,2,2);
  plot(Y_eval2(i,1:2:end),Y_eval2(i,2:2:end),'bo-',fX_eval2(i,1:2:end),fX_eval2(i,2:2:end),'ro-');   
  title(['rbf: i=' num2str(i) '; e2=' num2str(mean(tonguermse(fX_eval2(i,:),Y_eval2(i,:))),'%1.3fmm')]); 
  axis([min(min(Y_eval2(:,1:2:end))) max(max(Y_eval2(:,1:2:end))) min(min(Y_eval2(:,2:2:end))) max(max(Y_eval2(:,2:2:end)))]); 
  set(gca,'DataAspectRatio',[1 1 1]); axis ij; drawnow; 
  pause(0.5)
end


