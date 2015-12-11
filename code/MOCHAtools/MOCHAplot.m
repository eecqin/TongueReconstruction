% MOCHAplot(E,A,L[,I,axis_xy,pal,pha,vt,mov,mr])
%
% Given an utterance in the MOCHA-TIMIT dataset:
% - Right panel of the figure: selectively display acoustical features,
%   phonetic labels and EMA trajectories.
% - Left panel of the figure: optionally plot the animated vocal tract, and
%   a the moving bar in the right subplots.
%
% Needs the following external functions:
% - Several functions from our speech_analysis toolbox (to compute acoustics).
% - Malcolm Slaney's MakeQTMovie.m (to create the QuickTime movie).
% 
% In:
%  E: 1x1 structure, the EMA trajectories and their sampling rate.
%  A: 1x1 structure, the acoustic wave, sampling rate and nbits information.
%  L: ?x1 structure, the phonetical labels information.
%  I: indices for desired plots, defined as:
%      1: 'LIx' EMA channel for lower incisor X
%      2: 'LIy' EMA channel for lower incisor Y
%      3: 'ULx' EMA channel for upper lip X
%      4: 'ULy' EMA channel for upper lip Y
%      5: 'LLx' EMA channel for lower lip X
%      6: 'LLy' EMA channel for lower lip Y
%      7: 'TTx' EMA channel for tongue tip X
%      8: 'TTy' EMA channel for tongue tip Y
%      9: 'TBx' EMA channel for tongue body X
%     10: 'TBy' EMA channel for tongue body Y
%     11: 'TDx' EMA channel for tongue dorsum X
%     12: 'TDy' EMA channel for tongue dorsum Y
%     13: 'Vx'  EMA channel for velum X
%     14: 'Vy'  EMA channel for velum Y
%     15: 'UIx' EMA channel for upper incisor X (reference coil)
%     16: 'UIy' EMA channel for upper incisor Y (reference coil)
%     17: 'BNx' EMA channel for bridge nose X (reference coil)
%     18: 'BNy' EMA channel for bridge nose Y (reference coil)
%      0: vocal tract animation and moving bar
%     -1: acoustic wave and phonetic label
%     -2: spectrogram and formants
%     -3: energy contour
%     -4: pitch contour
%     For example, I = [0 -1 -2 1 2] plots: acoustic wave and phonetic label, 
%     spectrogram and formants, vocal tract animation and moving bar, and
%     the {'LIx','LIy'} EMA channels.
%  axis_xy: 1x4 vector of display ranges for X and Y axes, i.e.,
%     [min_x max_x min_y max_y].
%  pal: ?x2 matrix of palate outline. The 1st and 2nd columns contain X and Y 
%     coordinates of samples of the palate outline, respectively. 
%  pha: ?x2 matrix of pharyngeal outline. 
%  vt: ?x2 matrix of vocal tract outline. 
%  mov: filename of QuickTime movie. Default: don't create movie.
%  mr: moving rate for the process bar (e.g. mr=2 displays the process bar
%     every 2 articulatory frames). Useful to create smaller files. Default: 1.
% 
% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2008 by Chao Qin and Miguel A. Carreira-Perpinan

function MOCHAplot(E,A,L,I,axis_xy,pal,pha,vt,mov,mr)

% --------------- Argument defaults ---------------- %
if ~exist('I','var') | isempty(I) I = [-1 1:14]; end;
if ~exist('axis_xy','var') | isempty(axis_xy) axis_xy = []; end;
if ~exist('pal','var') | isempty(pal) pal = []; end;
if ~exist('pha','var') | isempty(pha) pha = []; end;
if ~exist('vt','var') | isempty(vt) vt = []; end;
if ~exist('mov','var') | isempty(mov) mov = []; end;
if ~exist('mr','var') | isempty(mr) mr = 1; end;
% --------------- Argument defaults ---------------- %

% Check if inputs are empty
if isempty(E) | isempty(A) | isempty(L)
  error('Required inputs are empty!');
end

% Plot configurations
left_offset = 0.1;	% Space to the left (for the Y axis label)
right_offset = 0.02;	% Space to the right (for the Y axis label)
sep = 0.015;		% Separation between subplots
top_offset = 0.03;	% Space to the top
bottom_offset = 0.03;	% Space to the bottom
nSubplot = length(I(find(I~=0)));
height = (1-top_offset-bottom_offset)/nSubplot;	% Height of individual graphs

% Acoustic configurations
wintime = 0.025;
hoptime = 0.01;
nCep = 13;
ceplifter = 22;
cepN = 1;
preemph = 0.97;
window = @hamming;
minfreq = 0;
maxfreq = 8000;
nChan = 26;
nFFT = 512;
p = 20;			% for 16KHz

% EMA
ema = E.dat; 
fs_ema = E.fs;

% Acoustic wave
y = A.dat; 
fs_wav = A.fs;
nbits = A.nbits;
nSample = size(y,1);
win_size = round(wintime*fs_wav);	% in samples
frame_shift = round(hoptime*fs_wav);	% in samples
nOverlap = win_size - frame_shift;
nFrame = floor((nSample-nOverlap)/frame_shift);
time_axis = (1:nSample)/fs_wav;

if size(ema,1)/fs_ema<=nSample/fs_wav 
  max_t = size(ema,1)/fs_ema; 
else
  max_t = nSample/fs_wav;
end

k = 1;
% Plot channels of interest
if ismember(0,I)
  % plot vocal tract animation
  set(gcf,'Position',[0 0 1200 600],'PaperType','A4','Color',[1 1 1],...
    'PaperUnits','centimeters','PaperPosition',[0 0 16 16*sqrt(2)]);
  % set the width of left and right panels
  leftPanel_width = 0.4; rightPanel_width = 0.45;   
else
  set(gcf,'Position',[0 0 800/sqrt(2) 600],'PaperType','A4','Color',[1 1 1],...
    'PaperUnits','centimeters','PaperPosition',[0 0 16 16*sqrt(2)]);
  % set the width of left and right panels
  leftPanel_width = 0; rightPanel_width = 1-left_offset-right_offset; 
  left_offset = left_offset-sep-0.025;
end

yAxisLim = [];
if ismember(-1,I)
  subplot('position',[left_offset+leftPanel_width+sep+0.025 ...
    (nSubplot-k)*height+bottom_offset rightPanel_width height-sep]);
  % plot acoustic wave
  plot(time_axis,y*2^(nbits-1),'k'); hold on;
  % plot phonetic label
  for l=1:size(L,2)
    phone = L(l).phone; 
    ts = L(l).startt;
    te = L(l).endt;
    text(str2num(ts)+(str2num(te)-str2num(ts))/2,max(y*2^(nbits-1))*2-6000,...
         phone,'HorizontalAlignment','center','VerticalAlignment','top',...
         'FontSize',9);
    if ~strcmp(ts,'0')
      plot([str2num(ts) str2num(ts)],...
           [min(y*2^(nbits-1)) max(y*2^(nbits-1))*2],'b--');
    end
  end
  axis([0 max_t min(y*2^(nbits-1)) max(y*2^(nbits-1))*2]);
  set(gca,'XAxisLocation','top','Box','on','FontSize',10);
  yAxisLim = [yAxisLim; get(gca,'YLim')];
  hold off;  
  k = k + 1;
end

if ismember(-2,I)
  % Narrowband spectrogram
  pspectrumN = specgramN(y);
  freq_axis = (0:nFFT/2)/nFFT*fs_wav;
  subplot('position',[left_offset+leftPanel_width+sep+0.025 ...
    (nSubplot-k)*height+bottom_offset rightPanel_width height-sep]);
  image(time_axis,freq_axis/1000,pspectrumN.*1000); hold on;
  colormap(1-gray);
  axis xy;  
  % Formant tracking
  [formant,bandwidth] = formant_track(y,fs_wav,wintime,hoptime);
  plot((1:nFrame)'*hoptime,formant(:,1)/1000,'r');
  plot((1:nFrame)'*hoptime,formant(:,2)/1000,'g');
  plot((1:nFrame)'*hoptime,formant(:,3)/1000,'b');
  ylabel('kHz','FontSize',13,'Position',[-0.1 mean(get(gca,'YLim'))]);
  xlim([0 max_t]); 
  if k==1
    set(gca,'XAxisLocation','top','Box','on','FontSize',10);
  else 
    set(gca,'XTickLabel',[],'Box','on');
  end
  yAxisLim = [yAxisLim; get(gca,'YLim')];
  hold off;  
  k = k + 1;
end

if ismember(-3,I)
  % Energy contour
  e = energycontour(y*2^(nbits-1),fs_wav,wintime,hoptime);
  subplot('position',[left_offset+leftPanel_width+sep+0.025 ...
    (nSubplot-k)*height+bottom_offset rightPanel_width height-sep]);
  plot((1:nFrame)'*hoptime,e,'k-');  
  xlim([0 max_t]);
  ylabel('dB','FontSize',13,'Position',[-0.1 mean(get(gca,'YLim'))]);
  if k==1
    set(gca,'XAxisLocation','top','Box','on','FontSize',10);
  else
    set(gca,'XTickLabel',[],'Box','on');
  end
  yAxisLim = [yAxisLim; get(gca,'YLim')];
  k = k + 1;
end

if ismember(-4,I)
  % F0 by SIFT
  f0 = siftF0(y,fs_wav,wintime,hoptime);
  subplot('position',[left_offset+leftPanel_width+sep+0.025 ...
    (nSubplot-k)*height+bottom_offset rightPanel_width height-sep]);
  plot((1:nFrame)'*hoptime,f0,'r.');  
  axis([0 max_t 0 500]);
  ylabel('Hz','FontSize',13,'Position',[-0.1 mean(get(gca,'YLim'))]);
  if k==1
    set(gca,'XAxisLocation','top','Box','on','FontSize',10);
  else
    set(gca,'XTickLabel',[],'Box','on');
  end
  yAxisLim = [yAxisLim; get(gca,'YLim')];
  k = k + 1;
end
  
% Plot EMA trajectories
names = {'LIx','LIy','ULx','ULy','LLx','LLy',...
         'TTx','TTy','TBx','TBy','TDx','TDy',...
         'Vx','Vy',...
         'UIx','UIy','BNx','BNy'};
for i=I(find(I>0))
  subplot('position',[left_offset+leftPanel_width+sep+0.025 ...
    (nSubplot-k)*height+bottom_offset ...
    rightPanel_width height-sep]);
  plot((1:size(ema,1))'/E.fs,ema(:,i),'b'); hold on;
  xlim([0 max_t]);
  ylabel(names(i),'FontSize',13,'Position',[-0.1 mean(get(gca,'YLim'))]);
  if k==1
    set(gca,'XAxisLocation','top','Box','on','FontSize',10);
  else
    set(gca,'XTickLabel',[],'Box','on');
  end
  yAxisLim = [yAxisLim; get(gca,'YLim')];  
  k = k + 1;
end

if ismember(0,I)
  g = subplot('position',[left_offset-0.025 bottom_offset ...
        leftPanel_width 1-top_offset-bottom_offset-0.05]); 
  text(-30,30,['Time: 0 sec'],'FontSize',14); 
  h = axes('position',[left_offset+leftPanel_width+sep+0.025 ...
      bottom_offset rightPanel_width nSubplot*height-sep]);  
  plot([0 0],[-1 1],'r-');
  xlim([0 max_t]);
  axis off
  
  % Clean up memory
  clear E A L 
  
  % Create the QuickTime movie
  if ~isempty(mov)
    MakeQTMovie('start',mov);
    MakeQTMovie('quality',1);
  end
  
  for i=1:mr:size(ema,1)
    % plot vocal tract animation on the left panel
    delete(g);
    g = subplot('position',[left_offset-0.025 bottom_offset ...
      leftPanel_width 1-top_offset-bottom_offset-0.05]); hold on;
    title('MOCHA-TIMIT: Midsagittal Vocal Tract Animation','FontSize',14);
    text(50,20,['Time: ' num2str(i/fs_ema,'%1.3f') ' sec'],'FontSize',14); 
    MOCHAframe(ema(i,:),[],axis_xy,pal,pha);
    % plot process bar on the right panel
    delete(h);
    h = axes('position',[left_offset+leftPanel_width+sep+0.025 ...
      bottom_offset rightPanel_width nSubplot*height-sep]);  
    plot([(i-1)/fs_ema (i-1)/fs_ema],[-1 1],'r-');
    xlim([0 max_t]);
    axis off
    % Add frame to the video
    if ~isempty(mov) MakeQTMovie('addfigure'); end
  end
  % Complete the movie by adding the audio
  if ~isempty(mov)
    MakeQTMovie('framerate',ceil(size(ema,1)/mr/max_t));
    MakeQTMovie('addsound',y,fs_wav);
    MakeQTMovie('finish');
  end
  
end

