% MOCHAtraj(e[,I,axis_xy]) 
%
% Displays temporal trajectories for selected articulators
% 
% In:
%  e: 1x18 vector of original EMA data, i.e.
%     [LIx,LIy,ULx,ULy,LLx,LLy,TTx,TTy,TBx,TBy,TDx,TDy,Vx,Vy,UIx,UIy,BNx,BNy].
%  I: indices for EMA channels. Default: [1:1:7].
%  axis_xy: axis for time frame (x) and articulatory position (y)

% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2009 by Chao Qin and Miguel A. Carreira-Perpinan

function MOCHAtraj(e,I,axis_xy)

% --------------- Argument defaults ---------------- %
if ~exist('I','var') | isempty(I) I = [1:1:7]; end;
if ~exist('axis_xy','var') | isempty(axis_xy) axis_xy = [0 size(e,1) -40 80]; end;
% --------------- Argument defaults ---------------- %

% Plot temporal trajectories for selected articulators
II = [2*setdiff(I,[8 9])'-1 2*setdiff(I,[8 9])']; II = II';
e1 = e(:,II(:));
II = [2*intersect(I,[8 9])'-1 2*intersect(I,[8 9])']; II = II';
e2 = e(:,II(:));

names = {'LI','UL','LL','TT','TB','TD','V','UI','BN'};
legendNameX=[]; legendNameY=[];
for i=setdiff(I,[8 9])
  if i==7
    legendNameX = [legendNameX; names{i} 'x '];
    legendNameY = [legendNameY; names{i} 'y '];  
  else
    legendNameX = [legendNameX; names{i} 'x'];
    legendNameY = [legendNameY; names{i} 'y'];  
  end
end
legendName = [legendNameX; legendNameY];

plot(e1(:,[1:2:size(e1,2)-1]),'-','LineWidth',1); hold on;
plot(e1(:,[2:2:size(e1,2)]),'--','LineWidth',1); 
if ismember(8,I)
  plot(e2(:,1),'k-','LineWidth',1); 
  plot(e2(:,2),'k--','LineWidth',1); 
  legendName = [legendName; 'UIx'; 'UIy'];
end
if ismember(9,I)
  plot(e2(:,3),'k-','LineWidth',1); 
  plot(e2(:,4),'k--','LineWidth',1); 
  legendName = [legendName; 'BNx'; 'BNy'];
end
set(gca,'FontSize',12)
xlabel('Time frame','FontSize',12)
ylabel('mm','FontSize',12)
axis(axis_xy);
legend(legendName,'Location','BestOutSide');
drawnow; hold off;
