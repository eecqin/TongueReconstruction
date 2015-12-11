format compact


clc; clear all; close all
K = 3; P = 20; N = 100; M = 20;
x = linspace(0,2,P); y = sin((1+randn(N,1)/5)*x+repmat(randn(N,1)/10,1,P));
x = repmat(x,N,1)+randn(N,P)/40;
Y1 = [];for i=1:N tmp = [x(i,:);y(i,:)]; Y1 = [Y1;tmp(:)']; end; Y1 = 70*Y1;
X1 = Y1(:,[5 6 11 12 31 32]); % K = 3;

% 2D-wise transformation
iA0 = 0.5*[cos(0.7) -sin(0.7);sin(0.7) cos(0.7)]; ib0 = [-100;50];
A0 = inv(iA0); b0 = -inv(iA0)*ib0; % This is the ground truth
Y2 = (kron(eye(P),iA0)*Y1')'+repmat(repmat(ib0',1,P),N,1);
X2 = Y2(:,[5 6 11 12 31 32]);

l = 0; [f,fX,E] = rbftrain(X1,Y1,M,[],l);	% Fit and prediction error
l = 0; [fl,flX,El] = linftrain(X1,Y1,l);	% Fit and prediction error
A = eye(2); b = [0;0];	% Initial A, b
Ab = [repmat(A(:),K,1); repmat(b(:),K,1); repmat(A(:),P,1); repmat(b(:),P,1)]'; % for F
Ab2 = [repmat(A(:),K,1); repmat(b(:),K,1); repmat(A(:),P,1); repmat(b(:),P,1)]'; % for E
Ab0 = [repmat(A0(:),K,1); repmat(b0(:),K,1); repmat(A0(:),P,1); repmat(b0(:),P,1)]';

% test gradients of rbf case
l = 0;
paramf = {f,X1,Y1,l};
% obj fcn F
[ff,gg,HH] = Ferr_rbf2(Ab,f,X1,Y1);
g = numgradhess(@Ferr_rbf2,paramf,Ab); % Numerical gradient & Hessian
% obj fcn E
[ff2, gg2] = Ferr2_rbf2(Ab,f,X1,Y1,l);
g2 = numgradhess(@Ferr2_rbf2,paramf,Ab); % Numerical gradient & Hessian
% obj fcn E but only (Ax,bx) are opt variables
% $$$ [ff2a, gg2a] = Ferr2_rbf2a([repmat(A(:),K,1); repmat(b(:),K,1)]',f,X1,Y1);
% $$$ g2a = numgradhess(@Ferr2_rbf2a,{f,X1,Y1},[repmat(A(:),K,1); repmat(b(:),K,1)]'); % Numerical gradient & Hessian


% test gradients of linf case
paramfl = {fl, X1, Y1};
% obj fcn F
[ffl, ggl] = Ferr_linf2(Ab,fl,X1,Y1);
gl = numgradhess(@Ferr_linf2,paramfl,Ab); % Numerical gradient & Hessian
% obj fcn E
[ffl2, ggl2] = Ferr2_linf2(Ab,fl,X1,Y1);
gl2 = numgradhess(@Ferr2_linf2,paramfl,Ab); % Numerical gradient & Hessian
% Ab is optimal 
tol = 1e-2;  % set this carefully
Ab_2 = rbfadapt2(X1,Y1,f,Ab',l,tol);   % optimize over obj fcn F: Ab2 should be close to Ab
Ab2_2 = rbfadapt22(X1,Y1,f,Ab2',l,tol); % optimize over obj fcn E: Ab22 should be close to Ab




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% 2D-wise transformation
A = eye(2); b = [0;0];	% Initial A, b
Ab = [repmat(A(:),K,1); zeros(2*K,1); repmat(A(:),P,1); zeros(2*P,1)]';
Ab2 = [repmat(A(:),K,1); zeros(2*K,1); repmat(A(:),P,1); zeros(2*P,1)]';
A0 = inv(params.adapt.iA0); b0 = -inv(params.adapt.iA0)*params.adapt.ib0; % This is the ground truth
Ab0 = [repmat(A0(:),K,1); repmat(b0(:),K,1); repmat(A0(:),P,1); repmat(b0(:),P,1)]';
Ab02 = [repmat(A0(:),K,1); repmat(b0(:),K,1); repmat(params.adapt.iA0(:),P,1); repmat(params.adapt.ib0(:),P,1)]';

I = [2 9 19]'; J = [1:24]';
II = [2*I-1 2*I]'; JJ = [2*J-1 2*J]';
yTR = V(iTR,JJ(:)); xTR = yTR(:,II(:)); yTE = V(iTE,JJ(:)); xTE = yTE(:,II(:));
Y_adpt = yTE(1:50,:); X_adpt = Y_adpt(:,II(:)); 
Y_adpt2 = (kron(eye(P),params.adapt.iA0)*Y_adpt')'+repmat(repmat(params.adapt.ib0',1,P),size(Y_adpt,1),1); X_adpt2 = Y_adpt2(:,II(:));

% --- MODEL
fl = linftrain(xTR,yTR);
load(['rbfcentres' num2str(params.rbf.M) '.mat'],'mu','-mat'); 
muI = mu(:,II(:)); % kmeans
f = rbftrain(xTR,yTR,muI,params.rbf.sigma,params.rbf.l);	 % Fit rbf  
% --- CORE
tol = 1e-3;
tic; Ab_2 = rbfadapt2(xTR(1:100,:),yTR(1:100,:),f,Ab',0,tol); toc;		% Adaptation on F
tic; Ab2_2 = rbfadapt22(xTR(1:100,:),yTR(1:100,:),f,Ab2',0,3.5); toc;		% Adaptation on E
tic; Ab0_2 = rbfadapt2(X_adpt2,Y_adpt2,f,Ab0',0,tol); toc;		% Adaptation on F
tic; Ab02_2 = rbfadapt22(X_adpt2(1:10,:),Y_adpt2(1:10,:),f,Ab02',0,0.05); toc;		% Adaptation on E

[Y,fgX,gX] = adapted_f2(xTR(1:100,:),f,Ab);  % fcn F
[Y_2,fgX_2,gX_2] = adapted_f2(xTR(1:100,:),f,Ab_2); % fcn F


[Y2,fgX2,gX2] = adapted_f22(xTR(1:100,:),f,Ab2);  % fcn E
[Y2_2,fgX2_2,gX2_2] = adapted_f22(xTR(1:10,:),f,Ab2_2); % fcn E
[Y2,fgX2,gX2] = adapted_f22(X_adpt2(1:10,:),f,Ab02);  % fcn E
[Y2_2,fgX2_2,gX2_2] = adapted_f22(X_adpt2(1:10,:),f,Ab02_2); % fcn E

figure(1);clf;
subplot(121); plot(Ab,'b-');hold on;plot(Ab_2,'r-');title('BFGS on fcn F');
subplot(122); plot(Ab2,'b-');hold on;plot(Ab2_2,'r-');title('BFGS on fcn E');
figure(2);
for i=1:100
  set(0,'CurrentFigure',2);clf;
  subplot(121);
  plot(yTR(i,1:2:end-1),yTR(i,2:2:end),'bo');hold on;
  plot(Y(i,1:2:end-1),Y(i,2:2:end),'ro');
  plot(Y_2(i,1:2:end-1),Y_2(i,2:2:end),'ko');
  axis ij; set(gca,'DataAspectRatio',[1 1 1]); 
  axis([40 130 40 100]);
  title(['i=' num2str(i)  '  Ab: ' num2str(norm(Y(i,:)-yTR(i,:)),3) '   Ab\_2: ' num2str(norm(Y_2(i,:)-yTR(i,:)),3)]);
  subplot(122);
  plot(yTR(i,1:2:end-1),yTR(i,2:2:end),'bo');hold on;
  plot(Y2(i,1:2:end-1),Y2(i,2:2:end),'ro');
  plot(Y2_2(i,1:2:end-1),Y2_2(i,2:2:end),'ko');
  axis ij; set(gca,'DataAspectRatio',[1 1 1]); 
  axis([40 130 40 100]);
  title(['i=' num2str(i)  '  Ab2: ' num2str(norm(Y2(i,:)-yTR(i,:)),3) '   Ab2\_2: ' num2str(norm(Y2_2(i,:)-yTR(i,:)),3)]);
  drawnow; pause
end

figure(5);
for i=1:10
  set(0,'CurrentFigure',5);clf;
  plot(Y_adpt2(i,1:2:end-1),Y_adpt2(i,2:2:end),'bo');hold on;
  plot(Y_2(i,1:2:end-1),Y_2(i,2:2:end),'ro');
  plot(Y2_2(i,1:2:end-1),Y2_2(i,2:2:end),'ko');
  axis ij; set(gca,'DataAspectRatio',[1 1 1]); 
  drawnow; pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NN = 10; tol = 1;
tic; [Al,bl] = linfadapt(X_adpt2(1:NN,:),Y_adpt2(1:NN,:),fl,0); toc;        % fcn F for nonadaptive case by linf
tic; [A,b] = rbfadapt(X_adpt2(1:NN,:),Y_adpt2(1:NN,:),f,A0,b0,0,0.05); toc; % fcn F for nonadaptive case rbf
iA = inv(A); ib = -inv(A)*b;
Ab0 = [repmat(A(:),K,1); repmat(b(:),K,1); repmat(iA(:),P,1); repmat(ib(:),P,1)]';
% fcn E for adaptive case by rbf+BFGS with initialization from rbfadapt
tic; Ab02_2 = rbfadapt22(X_adpt2(1:NN,:),Y_adpt2(1:NN,:),f,Ab0',0,20); toc;
% fcn E for adaptive case by rbf+AlterOpt with initialization from rbfadapt
tic; Ab02_2a = rbfadapt22a(X_adpt2(1:NN,:),Y_adpt2(1:NN,:),f,Ab0',1e-1,5); toc;

[Y,fgX,gX] = adapted_f(X_adpt2(1:NN,:),f,A,b); % fcn F for nonadaptive case
[Y2_2,fgX2_2,gX2_2] = adapted_f22(X_adpt2(1:NN,:),f,Ab02_2); % fcn E for adaptive case by rbf+BFGS
[Y2_2a,fgX2_2a,gX2_2a] = adapted_f22(X_adpt2(1:NN,:),f,Ab02_2a); % fcn E for adaptive case by rbf+AlterOpt

figure(123);clf;
subplot(311);
plot(Y_adpt2(1,1:2:end-1),Y_adpt2(1,2:2:end),'bo-');hold on;
plot(Y(1,1:2:end-1),Y(1,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
title(num2str(norm(Y-Y_adpt2(1:NN,:))));
subplot(312);
plot(Y_adpt2(1,1:2:end-1),Y_adpt2(1,2:2:end),'bo-');hold on;
plot(Y2_2(1,1:2:end-1),Y2_2(1,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
title(num2str(norm(Y2_2-Y_adpt2(1:NN,:))));
subplot(313);
plot(Y_adpt2(1,1:2:end-1),Y_adpt2(1,2:2:end),'bo-');hold on;
plot(Y2_2a(1,1:2:end-1),Y2_2a(1,2:2:end),'ro-');set(gca,'DataAspectRatio',[1 1 1]);axis ij;
title(num2str(norm(Y2_2a-Y_adpt2(1:NN,:))));


figure(124);clf;
plot(Ab0,'b-');hold on;plot(Ab02_2,'r-');