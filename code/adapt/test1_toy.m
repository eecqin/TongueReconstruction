format compact


K = 3; P = 20; N = 100; M = 10;
x = linspace(0,2,P); y = sin((1+randn(N,1)/5)*x+repmat(randn(N,1)/10,1,P));
x = repmat(x,N,1)+randn(N,P)/40;
Y1 = [];for i=1:N tmp = [x(i,:);y(i,:)]; Y1 = [Y1;tmp(:)']; end; Y1 = 70*Y1;
X1 = Y1(:,[5 6 15 16 35 36]); % K = 3;

% $$$ figure(1);clf;
% $$$ subplot(121); plot(Y1(:,1:2:end-1),Y1(:,2:2:end),'bo-');hold on;
% $$$ set(gca,'DataAspectRatio',[1 1 1]);
% $$$ % Centre and scale dataset
% $$$ XX = Y1(:,[1:2:end-1]); YY = Y1(:,[2:2:end]);
% $$$ [UU,VV,mm] = pca([XX(:) YY(:)],2);
% $$$ % $$$ Y1 = (Y1-repmat(mm,size(Y1,1),P))./max(sqrt(VV'));
% $$$ Y1 = (Y1-repmat(mm,size(Y1,1),P))./repmat(sqrt(VV'),size(Y1,1),P);
% $$$ X1 = Y1(:,[5 6 15 16 35 36]);
% $$$ subplot(122); plot(Y1(:,1:2:end-1),Y1(:,2:2:end),'bo-');hold on;
% $$$ set(gca,'DataAspectRatio',[1 1 1]);

% 2D-wise transformation
iA0 = eye(2); ib0 = [0;0];
A0 = inv(iA0); b0 = -inv(iA0)*ib0; % This is the ground truth
Y2 = (kron(eye(P),iA0)*Y1')'+repmat(repmat(ib0',1,P),N,1);
X2 = Y2(:,[5 6 15 16 35 36]);

l = 0; [f,fX,E] = rbftrain(X1,Y1,M,[],l);	% Fit and prediction error
l = 0; [fl,flX,El] = linftrain(X1,Y1);	% Fit and prediction error
A = eye(2); b = [0;0];	% Initial A, b
Ab = [repmat(A(:),K,1); repmat(b(:),K,1); repmat(A(:),P,1); repmat(b(:),P,1)]'; % for F
Ab2 = [repmat(A(:),K,1); repmat(b(:),K,1); repmat(A(:),P,1); repmat(b(:),P,1)]'; % for E
Ab0 = [repmat(A0(:),K,1); repmat(b0(:),K,1); repmat(A0(:),P,1); repmat(b0(:),P,1)]';

% --- CORE
l = 0;
% fcn F for nonadaptive case by linf
tic; [Al,bl] = linfadapt(X2,Y2,fl,l); toc;
% fcn F for nonadaptive case rbf
tic; [A,b] = rbfadapt(X2,Y2,f,A0,b0,l,1e-3); toc; 

% initialization for adaptive cases
Abl = [repmat(Al(:),K,1); repmat(bl(:),K,1); repmat(Al(:),P,1); repmat(bl(:),P,1)]';
iAl = inv(Al); ibl = -inv(Al)*bl; Abl_2 = [repmat(Al(:),K,1); repmat(bl(:),K,1); repmat(iAl(:),P,1); repmat(ibl(:),P,1)]';
A = eye(2); b = [0;0];
Ab = [repmat(A(:),K,1); repmat(b(:),K,1); repmat(A(:),P,1); repmat(b(:),P,1)]';
iA = inv(A); ib = -inv(A)*b; Ab_2 = [repmat(A(:),K,1); repmat(b(:),K,1); repmat(iA(:),P,1); repmat(ib(:),P,1)]';
% end of initialization for adaptive cases

% $$$ %%%%%%%%% fcn F for adaptive case by linf
% $$$ tic; [Abl2,El2] = linfadapt2(X2,Y2,fl,l); toc;
% $$$ %%%%%%%%% fcn F for adaptive case by rbf+BFGS with initialization from rbfadapt
% $$$ tic; [Ab2,E2] = rbfadapt2(X2,Y2,f,Ab',l,1e-5); toc;

% fcn E for adaptive case by linf+BFGS with initialization from rbfadapt
tic; Abl2_2 = linfadapt22(X2,Y2,fl,Abl_2',l,1e-6); toc;
% fcn E for adaptive case by rbf+BFGS with initialization from rbfadapt
tic; Ab2_2 = rbfadapt22(X2,Y2,f,Ab_2',l,1e-2); toc;
% fcn E for adaptive case by linf+AlterOpt with initialization from rbfadapt
tic; Abl2_2a = rbfadapt22a(X2,Y2,fl,Abl_2',1e-6,500); toc;
% fcn E for adaptive case by rbf+AlterOpt with initialization from rbfadapt
tic; Ab2_2a = rbfadapt22a(X2,Y2,f,Ab_2',1e-3,500); toc;

[Ab_2(1:12)' Abl2_2(1:12)' Ab2_2(1:12)' Abl2_2a(1:12)' Ab2_2a(1:12)']

[Yl,fgXl,gXl] = adapted_f(X2,fl,Al,bl);              % fcn F for nonadaptive case by linf
[Y,fgX,gX] = adapted_f(X2,f,A,b);                    % fcn F for nonadaptive case by rbf
[Yl2_2,fgXl2_2,gXl2_2] = adapted_f22(X2,fl,Abl2_2);  % fcn E for adaptive case by linf+BFGS
[Y2_2,fgX2_2,gX2_2] = adapted_f22(X2,f,Ab2_2);       % fcn E for adaptive case by rbf+BFGS
[Yl2_2a,fgXl2_2a,gXl2_2a] = adapted_f22(X2,fl,Abl2_2a);  % fcn E for adaptive case by linf+AlterOpt
[Y2_2a,fgX2_2a,gX2_2a] = adapted_f22(X2,f,Ab2_2a);   % fcn E for adaptive case by rbf+AlterOpt

i = 1;
figure(123);clf;
subplot(321);
plot(Y2(i,1:2:end-1),Y2(i,2:2:end),'bo-');hold on;
plot(Yl(i,1:2:end-1),Yl(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
title(['linfadapt on fcn F for nonadaptive case: ' num2str(norm(Yl-Y2))]);
subplot(322);
plot(Y2(i,1:2:end-1),Y2(i,2:2:end),'bo-');hold on;
plot(Y(i,1:2:end-1),Y(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
title(['rbfadapt (BFBS) on fcn F for nonadaptive case: ' num2str(norm(Y-Y2))]);
subplot(323);
plot(Y2(i,1:2:end-1),Y2(i,2:2:end),'bo-');hold on;
plot(Yl2_2(i,1:2:end-1),Yl2_2(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
title(['linfadapt (BFBS) on fcn E for adaptive case: ' num2str(norm(Yl2_2-Y2))]);
subplot(324);
plot(Y2(i,1:2:end-1),Y2(i,2:2:end),'bo-');hold on;
plot(Y2_2(i,1:2:end-1),Y2_2(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
title(['rbfadapt (BFGS) on fcn E for adaptive case: ' num2str(norm(Y2_2-Y2))]);
subplot(325);
plot(Y2(i,1:2:end-1),Y2(i,2:2:end),'bo-');hold on;
plot(Yl2_2a(i,1:2:end-1),Yl2_2a(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
title(['linfadapt (AlterOpt) on fcn E for adaptive case: ' num2str(norm(Yl2_2a-Y2))]);
subplot(326);
plot(Y2(i,1:2:end-1),Y2(i,2:2:end),'bo-');hold on;
plot(Y2_2a(i,1:2:end-1),Y2_2a(i,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
title(['rbfadapt (AlterOpt) on fcn E for adaptive case: ' num2str(norm(Y2_2a-Y2))]);

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
cond_Abl = ; Ax_Abl = Abl(1:4*K); Cy_Abl = Abl(6*K+1:6*K+4*P);
cond_Abl2_2 = ; Ax_Abl2_2 = Abl2_2(1:4*K); Cy_Abl2_2 = Abl2_2(6*K+1:6*K+4*P);
cond_Abl2_2a = ; Ax_Abl2_2a = Abl2_2a(1:4*K); Cy_Abl2_2a = Abl2_2a(6*K+1:6*K+4*P);
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