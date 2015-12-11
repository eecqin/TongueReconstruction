load maaw0_adapt.mat
V = s1.V; iTR = s1.iTR; iTE = s1.iTE;
I = [2 9 19]'; J = [1:24]';
II = [2*I-1 2*I]'; JJ = [2*J-1 2*J]';
yTR = V(iTR,JJ(:)); xTR = yTR(:,II(:));
yTE = V(iTE,JJ(:)); xTE = yTE(:,II(:));

% Nonadaptive 2D-wise transformation
K = size(xTR,2)/2; P = size(yTR,2)/2;

blka = zeros (2*P, 2*P);
ba = zeros(2*P,1);
Ab = zeros(4*P,1);
iA0 = [0.35 -2.5; 0.5 0.35]; ib0 = [150; -50];

% for p=1:P
%    blka(2*p-1:2*p,2*p-1:2*p) = iA0 + 0.01*randn(2,2); % 
%    ba (2*p-1:2*p) = ib0 + 0.01*randn(2,1);
% end

teta1 = (15)*pi/180; teta2 = (10)*pi/180; teta3 = (-10)*pi/180; teta4 = (-5)*pi/180;
ib0 = [20;-20];
iA1 = [cos(teta1) -sin(teta1); sin(teta1) cos(teta1)]; ib1 = ib0 + [-5;0];
iA2 = [cos(teta2) -sin(teta2); sin(teta2) cos(teta2)]; ib2 = ib0 + [-15;7];
iA3 = [cos(teta3) -sin(teta3); sin(teta3) cos(teta3)]; ib3 = ib0 + [-35;35];
iA4 = [cos(teta4) -sin(teta4); sin(teta4) cos(teta4)]; ib4 = ib0 + [-30;30];
for p=1:6
   blka(2*p-1:2*p,2*p-1:2*p) = iA1; % 
   ba (2*p-1:2*p) = ib1;
   Ab (4*(p-1)+1:4*p) = iA1(:);
end
for p=7:12
   blka(2*p-1:2*p,2*p-1:2*p) = iA2; % 
   ba (2*p-1:2*p) = ib2;
   Ab (4*(p-1)+1:4*p) = iA2(:);
end
for p=13:18
   blka(2*p-1:2*p,2*p-1:2*p) = iA3; % 
   ba (2*p-1:2*p) = ib3;
   Ab (4*(p-1)+1:4*p) = iA3(:);
end
for p=19:24
   blka(2*p-1:2*p,2*p-1:2*p) = iA4; % 
   ba (2*p-1:2*p) = ib4;
   Ab (4*(p-1)+1:4*p) = iA4(:);
end

Ab = [Ab;ba];
Y22 = trans(yTR, Ab,0);

figure(123);clf;
for i=1:10:size(yTR,1)
    Y = yTR(i,:); X = xTR(i,:); K = size(X,2)/2; P = size(Y,2)/2;
    Y2 = (blka*Y')'+ ba'; X2 = Y2(:,II(:));
    clf; set(0,'CurrentFigure',123); 
    plot(Y(1:2:end-1),Y(2:2:end),'bo-');hold on;
    plot(Y2(1:2:end-1),Y2(2:2:end),'ro-');
    set(gca,'DataAspectRatio',[1 1 1]);axis ij;
    legend('old speaker','new speaker','Location','NorthEast');
    drawnow;pause;
end