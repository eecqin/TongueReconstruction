% XRMBtraj(X,I,findex,title_text,yAxisRange)
%
% 
% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2009 by Chao Qin

function XRMBtraj(X,I,findex,title_text,yAxisRange)

% --------------- Argument defaults ---------------- %
if ~exist('I','var') | isempty(I) I = [1:16]; end; 
if ~exist('findex','var') | isempty(findex) findex = [1 size(X,1)]; end; 
if ~exist('title_text','var') | isempty(title_text) title_text = []; end;
if ~exist('yAxisRange','var') | isempty(yAxisRange) yAxisRange = []; end;
% --------------- Argument defaults ---------------- %
  
% $$$ [N,D] = size(X.true);
[N,D] = size(X);
  
% Plot configurations
left_offset = 0.1;		% Space to the left (for the Y axis label)
right_offset = 0.02;            % Space to the right (for the Y axis label)
sep = 0.01;			% Separation between subplots
top_offset = 0.03;	        % Space to the top
bottom_offset = 0.03;           % Space to the bottom
nSubplot = length(I);
height = (1-top_offset-bottom_offset)/nSubplot;  % Height of individual graphs

leftPanel_width = 0; rightPanel_width = 1-left_offset-right_offset; 
left_offset = left_offset-sep-0.025;

ARTnames = {'ULx','ULy',...
            'LLx','LLy',...
            'T1x','T1y',...
            'T2x','T2y',...
            'T3x','T3y',...
            'T4x','T4y',...
            'MNIx','MNIy',...
            'MNMx','MNMy'};


% Plot articulatory channels
% $$$ set(gcf,'Position',[0 0 600 60*length(I)],'PaperType','A4','Color',[1 1 1],'PaperUnits','centimeters','PaperPosition',[0 0 20 2*length(I)]);
kk=1;
for i=I
  subplot('position',[left_offset+leftPanel_width+sep+0.025 ...
    (nSubplot-kk)*height+bottom_offset rightPanel_width height-sep]);
% $$$   plot((1:N)',X.true(:,i),'b-'); hold on;
% $$$   plot((1:N)',X.rec(:,i),'r-'); hold off;
  plot([findex(1):findex(end)]'/145.6,X([findex(1):findex(end)],i),'b-'); hold on;
  xlim([findex(1)/145.6 findex(end)/145.6]);
  ylim([floor(yAxisRange.min(i)) ceil(yAxisRange.max(i))]);
  yAxisMean = get(gca,'YLim'); yAxisMean = mean(yAxisMean);
  ylabel(ARTnames(i),'FontSize',12);%,'Position',[-100 yAxisMean 1]);
% $$$   set(gca,'XAxisLocation','top');
  if kk==1 title(title_text,'FontSize',12); end;
  if kk~=I(end) set(gca,'XTickLabel',[],'Box','on'); drawnow; end;
  kk = kk+1; 
end


