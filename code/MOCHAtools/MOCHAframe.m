% MOCHAframe(e[,I,axis_xy,pal,pha,vt]) 
%
% Displays selected articulators and an interpolated tongue contour in the
% midsaggital plane for a single frame. The tongue is interpolated by a
% cubic B-spline (Matlab function 'spline'); specifically, we consider a
% uniform grid of P locations along the X axis (with known Y values for 3
% tongue pellets) and pass it to the spline function. Note that interpolation
% could be replaced by prediction which could yield more realistic tongue
% shapes (Qin et al: "Predicting tongue shapes from a few landmark locations",
% Interspeech 2008).
% 
% In:
%  e: 1x18 vector of original EMA data, i.e.
%     [LIx,LIy,ULx,ULy,LLx,LLy,TTx,TTy,TBx,TBy,TDx,TDy,Vx,Vy,UIx,UIy,BNx,BNy].
%  I: indices for EMA channels. Default: [1:1:7].
%  axis_xy, pal, pha, vt: see MOCHAplot.

% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2009 by Chao Qin and Miguel A. Carreira-Perpinan

function MOCHAframe(e,I,axis_xy,pal,pha,vt)

% --------------- Argument defaults ---------------- %
if ~exist('I','var') | isempty(I) I = [1:1:7]; end;
if ~exist('axis_xy','var') | isempty(axis_xy) axis_xy = [-40 80 -50 30]; end;
if ~exist('pal','var') | isempty(pal) pal = []; end;
if ~exist('pha','var') | isempty(pha) pha = []; end;
if ~exist('vt','var') | isempty(vt) vt = []; end;
% --------------- Argument defaults ---------------- %

min_x = axis_xy(1); max_x = axis_xy(2);
min_y = axis_xy(3); max_y = axis_xy(4); 
e_x = e(1:2:end-1); e_y = e(2:2:end);

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

% Plot tongue pellets
if isequal(ismember([4 5 6],I),[1 1 1])
  if e_x(4)<=min_x || e_x(4)>=max_x || e_y(4)<=min_y || e_y(4)>=max_y
    % plot TB and TD
    plot(e_x([5 6]),e_y([5 6]),'ro-','MarkerSize',5,'MarkerFaceColor','r');
    hold on;
  elseif e_x(5)<=min_x || e_x(5)>=max_x || e_y(5)<=min_y || e_y(5)>=max_y
    % plot TT and TD
    plot(e_x([4 6]),e_y([4 6]),'ro-','MarkerSize',5,'MarkerFaceColor','r');
    hold on;
  elseif e_x(6)<=min_x || e_x(6)>=max_x || e_y(6)<=min_y || e_y(6)>=max_y
    % plot TB and TD
    plot(e_x([4 5]),e_y([4 5]),'ro-','MarkerSize',5,'MarkerFaceColor','r');
    hold on;
  else
    % spline interpolation of tongue pellets (TT,TB,TD)
    pp = spline(e_x(4:6),e_y(4:6));
    e_xx = linspace(min(e_x(4:6)),max(e_x(4:6)),100);
    e_yy = ppval(pp,e_xx);
    plot(e_xx,e_yy,'r-','LineWidth',3); hold on;
  end
end

% $$$ % temp
% $$$ title('MOCHA: Midsagittal Vocal Tract Animation for fsew0\_001','FontSize',14);hold on;
% $$$ text(-25,-5,['Mouth'],'FontSize',14); 
% $$$ plot([68 68],[-35 0],'k-','LineWidth',3);

% plot all articulators
plot(e_x(I),e_y(I),'ro','MarkerSize',8,'MarkerFaceColor','r');

% Add text labels for all articulators
names = {'LI','UL','LL','TT','TB','TD','V','UI','BN'};
for i=I
  text(e_x(i),e_y(i),names(i),'HorizontalAlignment','right','FontSize',14);
end
xlabel('X (mm)','FontSize',16); ylabel('Y (mm)','FontSize',16);
axis([min_x max_x min_y max_y]);
set(gca,'Box','on','DataAspectRatio',[1 1 1],'XTick',[min_x:20:max_x],...
  'YTick',[min_y:20:max_y],'FontSize',16);
drawnow; hold off;

