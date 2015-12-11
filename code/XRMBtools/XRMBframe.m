% XRMBframe(x[,axis_xy,pal,pha,vt]) 
%
% Displays 8 articulators and an interpolated tongue contour in the
% midsaggital plane for a single frame. The tongue is interpolated by a cubic
% B-spline (Matlab function 'spline'); specifically, we consider a uniform
% grid of P locations along the X axis (with known Y values for 3 tongue
% pellets) and pass it to the spline function. Note that interpolation could
% be replaced by prediction which could yield more realistic tongue shapes
% (Qin et al: "Predicting tongue shapes from a few landmark locations",
% Interspeech 2008).
%
% In:
%  x: 1x16 vector of original X-ray microbeam data, i.e., 
%     [ULx,ULy,LLx,LLy,T1x,T1y,T2x,T2y,T3x,T3y,T4x,T4y,MNIx,MNIy,MNMx,MNMy].
%  axis_xy, pal, pha, vt: see XRMBplot.

% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2008 by Chao Qin and Miguel A. Carreira-Perpinan

function XRMBframe(x,axis_xy,pal,pha)

% --------------- Argument defaults ---------------- %
if ~exist('axis_xy','var') | isempty(axis_xy) axis_xy = [-100 40 -40 40]; end;
if ~exist('pal','var') | isempty(pal) pal = []; end;
if ~exist('pha','var') | isempty(pha) pha = []; end;
if ~exist('vt','var') | isempty(vt) vt = []; end;
% --------------- Argument defaults ---------------- %

min_x = axis_xy(1); max_x = axis_xy(2);
min_y = axis_xy(3); max_y = axis_xy(4); 
x_x = x(1:2:15); x_y = x(2:2:16);

% Plot palate outline
if exist('pal','var') & ~isempty(pal)
  plot(pal(:,1),pal(:,2),'k-','LineWidth',3); hold on;
end

% Plot pharyngeal outline
if exist('pha','var') & ~isempty(pha)
  plot(pha(:,1),pha(:,2),'k-','LineWidth',3); hold on;
end

% Plot vocal tract outline
if exist('vt','var') & ~isempty(vt)
  plot(vt(:,1),vt(:,2),'k-','LineWidth',3); hold on;
end

% Plot articulators
if x_x(3)<=min_x || x_x(3)>=max_x || x_y(3)<=min_y || x_y(3)>=max_y
  % plot T2,T3,T4
  plot(x_x([4 5 6]),x_y([4 5 6]),'ro-','MarkerSize',5,'MarkerFaceColor','r');
  hold on;
elseif x_x(4)<=min_x || x_x(4)>=max_x || x_y(4)<=min_y || x_y(4)>=max_y
  % plot T1,T3,T4
  plot(x_x([3 5 6]),x_y([3 5 6]),'ro-','MarkerSize',5,'MarkerFaceColor','r');
  hold on;
elseif x_x(5)<=min_x || x_x(5)>=max_x || x_y(5)<=min_y || x_y(5)>=max_y
  % plot T1,T2,T4
  plot(x_x([3 4 6]),x_y([3 4 6]),'ro-','MarkerSize',5,'MarkerFaceColor','r');
  hold on;
elseif x_x(6)<=min_x || x_x(6)>=max_x || x_y(6)<=min_y || x_y(6)>=max_y
  % plot T1,T2,T3
  plot(x_x([3 4 5]),x_y([3 4 5]),'ro-','MarkerSize',5,'MarkerFaceColor','r');
  hold on;
else
  if sum(isnan([x_x(3:6) x_y(3:6)]))==0
    % spline interpolation of tongue pellets (TT,TB,TD)
    pp = spline(x_x(3:6),x_y(3:6));
    x_xx = linspace(min(x_x(3:6)),max(x_x(3:6)),100);
    x_yy = ppval(pp,x_xx);
    plot(x_xx,x_yy,'r-','LineWidth',3); hold on;
  end
  % plot T1,T2,T3,T4
  plot(x_x(3:6),x_y(3:6),'ro','MarkerSize',8,'MarkerFaceColor','r');
end
% plot UL, LL, MNI, and MNM
plot(x_x([1 2 7 8]),x_y([1 2 7 8]),'ro','MarkerSize',8,'MarkerFaceColor','r');

% Add text labels for all articulator
names = {'UL','LL','T1','T2','T3','T4','MNI','MNM'};
for i=1:8
  text(x_x(i),x_y(i),names(i),'HorizontalAlignment','right','FontSize',14);
end

xlabel('X (mm)','FontSize',16); ylabel('Y (mm)','FontSize',16);
axis([min_x max_x min_y max_y]);
set(gca,'Box','on','DataAspectRatio',[1 1 1],'XTick',[min_x:20:max_x],...
  'YTick',[min_y:20:max_y],'FontSize',16);
drawnow; hold off;

