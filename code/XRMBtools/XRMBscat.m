% XRMBscat(X[,axis_xy,pal,pha,vt])
%
% Scatterplot of X-ray microbeam data of 8 articulators in the midsaggital
% plane.
%
% In:
%  X: see XRMBtrace.
%  axis_xy, pal, pha, vt: see XRMBplot.

% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2008 by Chao Qin and Miguel A. Carreira-Perpinan

function XRMBscat(X,axis_xy,pal,pha,vt)

% --------------- Argument defaults ---------------- %
if ~exist('axis_xy','var') | isempty(axis_xy) axis_xy = [-100 40 -40 40]; end; 
if ~exist('pal','var') | isempty(pal) pal = []; end;
if ~exist('pha','var') | isempty(pha) pha = []; end;
if ~exist('vt','var') | isempty(vt) vt = []; end;
% --------------- Argument defaults ---------------- %

% Ranges of X-Y axes
min_x = axis_xy(1); max_x = axis_xy(2);
min_y = axis_xy(3); max_y = axis_xy(4); 

% Compute means and stdevs of each articulatory distribution
t = 0:2*pi/1000:2*pi;
xy = [cos(t)' sin(t)'];
xx = [0 1; -1 0; 0 -1; 1 0];
% $$$ art_cov = zeros(2,2,8);
% $$$ for i=1:8
% $$$   art_cov(:,:,i) = cov(X(:,[2*i-1 2*i]));
% $$$ end

% Plot X-ray microbeam data
colors = [0.5 0 1; 0 1 0; 1 0 1; 0.5 0.5 0; 1 1 0; 1 0 0; 0 1 1; 0.5 0.5 0.5];
pos = zeros(8,2);
for i=1:8
  plot(X(:,2*i-1),X(:,2*i),'.','Color',colors(i,:),'MarkerSize',1); hold on;
  pos(i,:) = mean(X(:,2*i-1:2*i)) - [5 10];
end

names = {'UL','LL','T1','T2','T3','T4','MNI','MNM'};
for i=1:8
  mu = mean(X(:,[2*i-1 2*i]));
% $$$   [V,D] = eig(art_cov(:,:,i));
  [V,D] = eig(cov(X(:,[2*i-1 2*i]))); 
  % $$$   [V,D] = pca(X(:,[2*i-1 2*i]),2); D = diag(D);
  xy_new = xy*(V*sqrt(D))' + repmat(mu,size(xy,1),1);
  yy = xx*(V*sqrt(D))' + repmat(mu,size(xx,1),1);
  plot([yy(1,1) yy(3,1)],[yy(1,2) yy(3,2)],'k-','LineWidth',1);
  plot([yy(2,1) yy(4,1)],[yy(2,2) yy(4,2)],'k-','LineWidth',1);
  plot(xy_new(:,1),xy_new(:,2),'k-','LineWidth',1);
% $$$   % Plot convex hull (optional)
% $$$   K = convhull(X(:,2*i-1),X(:,2*i));
% $$$   plot(X(K,2*i-1),X(K,2*i),'k-');
  % Add text labels for all articulator
  text(pos(i,1),pos(i,2),names(i),'FontSize',14); 
end
plot([-100 40],[0 0],'k-','LineWidth',1); 
plot([0 0],[-40 40],'k-','LineWidth',1); 

% Plot palate outline
if ~isempty(pal) plot(pal(:,1),pal(:,2),'k-','LineWidth',3); end

% Plot pharynx outline
if ~isempty(pha) plot(pha(:,1),pha(:,2),'k-','LineWidth',3); end

% Plot vocal tract outline
if ~isempty(vt) plot(vt(:,1),vt(:,2),'k-','LineWidth',3); end

xlabel('X (mm)','FontSize',16); ylabel('Y (mm)','FontSize',16);
axis([min_x max_x min_y max_y]);
set(gca,'Box','on','DataAspectRatio',[1 1 1],'XTick',[min_x:20:max_x],...
  'YTick',[min_y:20:max_y],'FontSize',16); hold off; drawnow; 

