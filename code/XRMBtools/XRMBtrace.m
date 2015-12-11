% XRMBtrace(X[,I,axis_xy,pal,pha,vt]) 
%
% Displays traces of selected articulators in the midsaggital plane 
% for a single utterance. 
% 
% In:
%  X: ?x16 matrix of X-ray microbeam data. Each row vector contains 
%     [ULx,ULy,LLx,LLy,T1x,T1y,T2x,T2y,T3x,T3y,T4x,T4y,MNIx,MNIy,MNMx,MNMy].
%  I: indices of desired articulatory channels. Default: [1 2 3 4 5 6 7 8],
%     i.e., {'UL','LL','T1','T2','T3','T4','MNI','MNM'}.
%  axis_xy, pal, pha, vt: see XRMBplot.

% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2008 by Chao Qin and Miguel A. Carreira-Perpinan

function XRMBtrace(X,I,axis_xy,pal,pha)

% --------------- Argument defaults ---------------- %
if ~exist('I','var') | isempty(I) I = [1:8]; end;
if ~exist('axis_xy','var') | isempty(axis_xy) axis_xy = [-100 40 -40 40]; end;
if ~exist('pal','var') | isempty(pal) pal = []; end;
if ~exist('pha','var') | isempty(pha) pha = []; end;
if ~exist('vt','var') | isempty(vt) vt = []; end;
% --------------- Argument defaults ---------------- %

% Ranges of X-Y axes
min_x = axis_xy(1); max_x = axis_xy(2);
min_y = axis_xy(3); max_y = axis_xy(4); 

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

% Plot traces of X-ray microbeam data
colors = [0.5 0 1; 0 1 0; 1 0 1; 0.5 0.5 0; 1 1 0; 1 0 0; 0 1 1; 0.5 0.5 0.5];
II = [2*I'-1 2*I']; II = II';
X2 = X(:,II(:));
for i=1:8
  plot(X2(:,2*i-1),X2(:,2*i),'.','Color',colors(i,:),'MarkerSize',2); hold on;
end

% Add text labels for all articulator
names = {'UL','LL','T1','T2','T3','T4','MNI','MNM'};
for i=I
  text(mean(X(:,2*i-1)),mean(X(:,2*i)),names(i),...
       'HorizontalAlignment','right','FontSize',14); 
end
axis([min_x max_x min_y max_y]);
xlabel('X (mm)','FontSize',16); ylabel('Y (mm)','FontSize',16);
set(gca,'Box','on','DataAspectRatio',[1 1 1],'XTick',[min_x:20:max_x],...
  'YTick',[min_y:20:max_y],'FontSize',16);
drawnow; hold off;

