function [hits,in0,in,x,y,xv0,yv0,xv,yv] = polyinterarea(fX,ui,li,v,pal,pha,db)
  
switch db
 case 'MOCHA'
  % first connect palate to velum with a spline
  pp = spline([pal(:,1); v(1)],[pal(:,2); v(2)]);
  e_xx = linspace(pal(end,1),v(1),10); e_yy = ppval(pp,e_xx);
  pal = [pal(1:end,:); [e_xx' e_yy']];
  
  [XX,YY] = meshgrid([-20:1:80],[-50:1:20]);
  xv0 = [fX([1:2:end-1]) fX(1)]'; yv0 = [fX([2:2:end]) fX(2)]';
  xv.pal = [pal(1,1); pal(:,1); pal(end,1); 30; pal(1,1)]; 
  yv.pal = [5; pal(:,2); 10; 10; 5];
  xv.pha = [pha(:,1); 80; 80; pha(1,1)]; 
  yv.pha = [pha(:,2); -40; -10; pha(1,2)]; 
  xv.ui = [ui(1)+2 ui(1) ui(1) ui(1) ui(1)+2 ui(1)+2]'; 
  yv.ui = [-4+ui(2) -4+ui(2) ui(2) 4+ui(2) 4+ui(2) -4+ui(2)]';
% $$$   xv.li = [li(1) li(1)+2 li(1)+2 li(1) li(1) li(1)]'; 
% $$$   yv.li = [li(2)-4 li(2)-4 li(2)+4 li(2)+4 li(2) li(2)-4]';
  xv.li = [li(1) li(1)-2 li(1)-2 li(1) li(1) li(1)]'; 
  yv.li = [li(2)-4 li(2)-4 li(2)+4 li(2)+4 li(2) li(2)-4]';
  [XX,YY] = meshgrid([-20:1:80],[-45:1:5]);
  x = XX(:); y = YY(:); 
  in0 = inpolygon(x,y,xv0,yv0); index0 = find(in0==1);
  in.pal = inpolygon(x,y,xv.pal,yv.pal); index.pal = find(in.pal==1);
  in.pha = inpolygon(x,y,xv.pha,yv.pha); index.pha = find(in.pha==1);
  in.ui = inpolygon(x,y,xv.ui,yv.ui); index.ui = find(in.ui==1);
  in.li = inpolygon(x,y,xv.li,yv.li); index.li = find(in.li==1);
 case 'XRMB'
  xv0 = [fX([1:2:end-1]) fX(1)]'; yv0 = [fX([2:2:end]) fX(2)]';
  xv.pal = [pal(1,1); pal(:,1); pal(end,1); -35; -4]; 
  yv.pal = [15; pal(:,2); 38; 35; 15];
  xv.pha = [pha(:,1); -100; -100; pha(1,1)]; 
  yv.pha = [pha(:,2); pha(2,2); pha(1,2); pha(1,2)];
  xv.ui = [ui(1)-2 ui(1) ui(1) ui(1) ui(1)-2 ui(1)-2]'; 
  yv.ui = [ui(2)+2 ui(2)+2 ui(2) ui(2)+10 ui(2)+10 ui(2)+2]';
  xv.ui = [ui(1)-2 ui(1) ui(1) ui(1) ui(1)-2 ui(1)-2]'; 
  yv.ui = [ui(2)-2 ui(2)-2 ui(2) ui(2)+8 ui(2)+8 ui(2)-2]';
  xv.li = [li(1) li(1)-2 li(1)-2 li(1) li(1) li(1)]'; 
  yv.li = [li(2)-4 li(2)-4 li(2)+4 li(2)+4 li(2) li(2)-4]';     
  [XX,YY] = meshgrid([-100:1:10],[-30:1:30]);
  x = XX(:); y = YY(:); 
  in0 = inpolygon(x,y,xv0,yv0); index0 = find(in0==1);
  in.pal = inpolygon(x,y,xv.pal,yv.pal); index.pal = find(in.pal==1);
  in.pha = inpolygon(x,y,xv.pha,yv.pha); index.pha = find(in.pha==1);
  in.ui = inpolygon(x,y,xv.ui,yv.ui); index.ui = find(in.ui==1);
  in.li = inpolygon(x,y,xv.li,yv.li); index.li = find(in.li==1);
end
hits.pal = intersect(index0,index.pal);
hits.pha = intersect(index0,index.pha);
hits.ui = intersect(index0,index.ui);
hits.li = intersect(index0,index.li);
