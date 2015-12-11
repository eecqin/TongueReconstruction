% MOCHAscat(E[,I,axis_xy,pal,pha,vt])
%
% Scatterplot of EMA data of 7 articulators in the midsaggital plane.
%
% In:
%  E: ?x18 matrix of EMA data for all utterances of one speaker. 
%  I, axis_xy, pal, pha, vt: see MOCHAtrace.

% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2008 by Chao Qin and Miguel A. Carreira-Perpinan

function MOCHAscat(E,I,axis_xy,pal,pha,vt)

% --------------- Argument defaults ---------------- %
if ~exist('I','var') | isempty(I) I = [1:1:7]; end;
if ~exist('axis_xy','var') | isempty(axis_xy) axis_xy = [-40 80 -50 30]; end;
if ~exist('pal','var') | isempty(pal) pal = []; end;
if ~exist('pha','var') | isempty(pha) pha = []; end;
if ~exist('vt','var') | isempty(vt) vt = []; end;
% --------------- Argument defaults ---------------- %

% Ranges of X-Y axes
min_x = axis_xy(1); max_x = axis_xy(2);
min_y = axis_xy(3); max_y = axis_xy(4); 

% Compute means and stdevs of each articulatory distribution
t = 0:2*pi/50:2*pi;
xy = [cos(t)' sin(t)'];
xx =[0 1; -1 0; 0 -1; 1 0];
art_cov = zeros(2,2,length(I));
for i=I
  art_cov(:,:,i) = cov(E(:,[2*i-1 2*i]));
end

% Plot palate outline
if exist('pal','var') & ~isempty(pal)
  plot(pal(:,1),pal(:,2),'k-','LineWidth',3); hold on;
end

% Plot pharynx outline
if exist('pha','var') & ~isempty(pha)
  plot(pha(:,1),pha(:,2),'k-','LineWidth',3); hold on;
end

% Plot vocal tract outline
if exist('vt','var') & ~isempty(vt)
  plot(vt(:,1),vt(:,2),'k-','LineWidth',3); hold on;
end

% Plot entire EMA data
II = [2*I'-1 2*I']; II = II';
E2 = E(:,II(:));
plot(E2(:,[1:2:end-1]),E2(:,[2:2:end]),'.','MarkerSize',0.5); hold on;
names = {'LI','UL','LL','TT','TB','TD','V','UI','BN'};
pos = [5 -35; -5 2; -8 -39; 18 -30; 32 -25; 45 -20; 62 -15; -7 -7; 10 50];
for i=I
  mu = mean(E(:,[2*i-1 2*i]));
  [V,D] = eig(art_cov(:,:,i));
  xy_new = xy*(V*sqrt(D))' + repmat(mu,size(xy,1),1);
  yy = xx*(V*sqrt(D))' + repmat(mu,size(xx,1),1);
  plot([yy(1,1) yy(3,1)],[yy(1,2) yy(3,2)],'k-','LineWidth',1);
  plot([yy(2,1) yy(4,1)],[yy(2,2) yy(4,2)],'k-','LineWidth',1);
  plot(xy_new(:,1),xy_new(:,2),'k-','LineWidth',1);
  % Add text labels for all articulator
  text(pos(i,1),pos(i,2),names(i),'FontSize',14);
end
plot([-40 100],[0 0],'k-','LineWidth',1); 
plot([0 0],[-70 70],'k-','LineWidth',1); 
xlabel('X (mm)','FontSize',16); ylabel('Y (mm)','FontSize',16);
axis([min_x max_x min_y max_y]);
set(gca,'Box','on','DataAspectRatio',[1 1 1],'XTick',[min_x:20:max_x],...
  'YTick',[min_y:20:max_y],'FontSize',16);
drawnow; hold off;

