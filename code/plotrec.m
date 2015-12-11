function plotrect2(fX,Y,L,pha,pal,db) % Simplified plot of plotrect.m

set(0,'CurrentFigure',gcf);clf;
h = text(0,0,'');
for i=1:size(fX,1)
  set(0,'CurrentFigure',1); clf; hold on;

  % plot palate and pharyngeal wall
  if exist('pha','var') & ~isempty(pha) plot(pha(:,1),pha(:,2),'k-','LineWidth',3); end
  if exist('pal','var') & ~isempty(pal) plot(pal(:,1),pal(:,2),'k-','LineWidth',3); end

  % plot original tongue pellets and the predicted tongue contour
  plot(Y(i,[1:2:end-1]),Y(i,[2:2:end]),'bo','MarkerSize',15,'LineWidth',2);
  plot(fX(i,[1:2:end-1]),fX(i,[2:2:end]),'ro-','MarkerSize',5,'MarkerFaceColor','r','LineWidth',2);
  % plot phonetic label
  switch db
   case 'XRMB'
    % spline interpolation of tongue pellets
    pp = spline(Y(i,1:2:end-1),Y(i,2:2:end));
    e_xx = linspace(min(fX(i,1:2:end-1)),max(fX(i,1:2:end-1))); e_yy = ppval(pp,e_xx);
    plot(e_xx,e_yy,'g-','LineWidth',2); 
    
    if exist('L','var') & ~isempty(L)
      delete(h); h = text(10,30.5,['Phone: /' L{i} '/'],'FontSize',15);
    end
    axis([-100 30 -40 40]); set(gca,'XTick',[-100:10:30],'YTick',[-40:10:40]);
   case 'MOCHA'
    % spline interpolation of tongue pellets
    pp = spline(Y(i,1:2:end-1),Y(i,2:2:end));
    e_xx = linspace(min(fX(i,1:2:end-1)),max(fX(i,1:2:end-1))); e_yy = ppval(pp,e_xx);
    plot(e_xx,e_yy,'g-','LineWidth',2);     					    
    if exist('L','var') & ~isempty(L)
      delete(h); h = text(58,13.5,['Phone: /' L{i} '/'],'FontSize',15);
    end
    axis([-20 80 -50 20]); set(gca,'XTick',[-20:10:80],'YTick',[-50:10:20]);
  end
  set(gca,'DataAspectRatio',[1 1 1],'Box','on','FontSize',16);
  xlabel('X (mm)'); ylabel('Y (mm)'); title(['Frame ' num2str(i)]); drawnow; pause
end

