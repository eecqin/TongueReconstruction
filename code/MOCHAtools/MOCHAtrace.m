% MOCHAtrace(E[,I,axis_xy,pal,pha,vt]) 
%
% Displays traces of selected articulators in the midsaggital plane 
% for a single utterance. 
% 
% In:
%  E: ?x18 matrix of EMA data for one utterance. Each row vector contains 
%     [LIx,LIy,ULx,ULy,LLx,LLy,TTx,TTy,TBx,TBy,TDx,TDy,Vx,Vy,UIx,UIy,BNx,BNy].
%  I: indices of desired articulatory channels. Default: [1 2 3 4 5 6 7], 
%     i.e. {'LI','UL','LL','TT','TB','TD','V'}.
%  axis_xy, pal, pha, vt: see MOCHAplot.

% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2008 by Chao Qin and Miguel A. Carreira-Perpinan

function MOCHAtrace(E,I,axis_xy,pal,pha)

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

% Plot traces of EMA data
II = [2*I'-1 2*I']; II = II';
E2 = E(:,II(:));
plot(E2(:,[1:2:end-1]),E2(:,[2:2:end]),'.','MarkerSize',1); hold on;

% Add text labels for all articulator
names = {'LI','UL','LL','TT','TB','TD','V','UI','BN'};
for i=I
  text(mean(E(:,2*i-1))-2,mean(E(:,2*i))+2,names(i),...
    'HorizontalAlignment','right','FontSize',14);
end
axis([min_x max_x min_y max_y]);
xlabel('X (mm)','FontSize',16); ylabel('Y (mm)','FontSize',16);
set(gca,'Box','on','DataAspectRatio',[1 1 1],'XTick',[min_x:20:max_x],...
  'YTick',[min_y:20:max_y],'FontSize',16);
drawnow; hold off;

