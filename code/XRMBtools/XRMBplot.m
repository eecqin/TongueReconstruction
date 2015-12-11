% XRMBplot(X,A,L[,I,axis_xy,pal,pha,vt,mov,mr])
%
% Given an utterance in the Wisconsin X-ray microbeam dataset:
% - Right panel of the figure: selectively display acoustical features,
%   phonetic labels and X-ray microbeam trajectories.
% - Left panel of the figure: optionally plot the animated vocal tract, and
%   a the moving bar in the right subplots.
%
% Needs the following external functions:
% - Several functions from our speech_analysis toolbox (to compute acoustics).
% - Malcolm Slaney's MakeQTMovie.m (to create the QuickTime movie).
%
% In:
%  X: 1x1 structure, the X-ray microbeam trajectories and their sampling rate,
%     and indices of good frames of X-ray microbeam recordings.
%  A: 1x1 structure, the acoustic wave, sampling rate and nbits information.
%  L: ?x1 structure, the phonetical labels information.
%  I: indices for desired plots, defined as:
%      1: 'ULx' X-ray channel
%      2: 'ULy' X-ray channel
%      3: 'LLx' X-ray channel
%      4: 'LLy' X-ray channel
%      5: 'T1x' X-ray channel
%      6: 'T1y' X-ray channel
%      7: 'T2x' X-ray channel
%      8: 'T2y' X-ray channel
%      9: 'T3x' X-ray channel
%     10: 'T3y' X-ray channel
%     11: 'T4x' X-ray channel
%     12: 'T4y' X-ray channel
%     13: 'MNIx' X-ray channel
%     14: 'MNIy' X-ray channel
%     15: 'MNMx' X-ray channel
%     16: 'MNMy' X-ray channel
%      0: vocal tract animation and moving bar
%     -1: acoustic wave and phonetic label
%     -2: spectrogram and formants
%     -3: energy contour
%     -4: pitch contour
%     For example, I = [0 -1 -2 1 2] plots: acoustic wave and phonetic label, 
%     spectrogram and formants, vocal tract animation and moving bar, and
%     the {'ULx','ULy'} X-ray channels.
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

function XRMBplot(X,A,L,I,axis_xy,pal,pha,vt,mov,mr)

% --------------- Argument defaults ---------------- %
if ~exist('I','var') | isempty(I) I = [-1 1:16]; end;
if ~exist('axis_xy','var') | isempty(axis_xy) axis_xy = []; end;
if ~exist('pal','var') | isempty(pal) pal = []; end;
if ~exist('pha','var') | isempty(pha) pha = []; end;
if ~exist('vt','var') | isempty(vt) vt = []; end;
if ~exist('mov','var') | isempty(mov) mov = []; end;
if ~exist('mr','var') | isempty(mr) mr = 1; end;
% --------------- Argument defaults ---------------- %

% Check if inputs are empty
if isempty(X) | isempty(A) error('Required inputs are empty!'); end

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
p = 20;           % for 16KHz

% X-ray microbeam
xray = X.dat; 
fs_xray = X.fs;

% Acoustic wave
y = A.dat; 
fs_wav = A.fs;
nbits = A.nbits;
nSample = size(y,1);
win_size = round(wintime*fs_wav);      % in samples
frame_shift = round(hoptime*fs_wav);   % in samples
nOverlap = win_size - frame_shift;
nFrame = floor((nSample-nOverlap)/frame_shift);
time_axis = (1:nSample)/fs_wav;

if size(xray,1)/fs_xray<=nSample/fs_wav 
  max_t = size(xray,1)/fs_xray; 
else
  max_t = nSample/fs_wav;
end

k = 1;
% Plot channels of interest
if ismember(0,I)
  % plot vocal tract animation
  set(gcf,'Position',[0 100 1200 600],'PaperType','A4','Color',[1 1 1],...
    'PaperUnits','centimeters','PaperPosition',[0 0 16 16*sqrt(2)]);
  % set the width of left and right panels
  leftPanel_width = 0.4; rightPanel_width = 0.45;   
else
  set(gcf,'Position',[0 100 800/sqrt(2) 600],'PaperType','A4','Color',[1 1 1],...
    'PaperUnits','centimeters','PaperPosition',[0 0 16 16*sqrt(2)]);
  % set the width of left and right panels
  leftPanel_width = 0; rightPanel_width = 1-left_offset-right_offset; 
  left_offset = left_offset-sep-0.025;
end

% Plot the right panel
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
    text(str2num(ts)+(str2num(te)-str2num(ts))/2,max(y*2^(nbits-1))+2500,...
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
  ylabel('kHz','FontSize',13,'Position',[-max_t/20 mean(get(gca,'YLim'))]);
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
  ylabel('dB','FontSize',13,'Position',[-max_t/20 mean(get(gca,'YLim'))]);
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
  % The code below is to compensate the mismatch of number of frames, i.e.,
  % the number of frames in f0 typically doesn't match the number of
  % articulatory frames. This is because in the SIFT algorithm there is a 
  % downsampling step which, combined with the stepsize (hoptime) in XRMB
  % (6.866 ms), results in the mismatch.
  % Thus, the function siftF0 here is used for demonstrating the plot. 
  % It is not a routine that implements the robust pitch extraction
  % algorithm. Users should use this function with caution and could replace 
  % it with other pitch extraction routine here. 
  if size(f0,1)<nFrame
    tmp = zeros(nFrame,1);
    tmp(1:size(f0,1)) = f0;
    f0 = tmp;
  end
  subplot('position',[left_offset+leftPanel_width+sep+0.025 ...
    (nSubplot-k)*height+bottom_offset rightPanel_width height-sep]);
  plot((1:nFrame)'*hoptime,f0,'r.');  
  axis([0 max_t 0 500]);
  ylabel('Hz','FontSize',13,'Position',[-max_t/20 mean(get(gca,'YLim'))]);
  if k==1
    set(gca,'XAxisLocation','top','Box','on','FontSize',10);
  else
    set(gca,'XTickLabel',[],'Box','on');
  end
  yAxisLim = [yAxisLim; get(gca,'YLim')];
  k = k + 1;
end
  
% Plot X-ray microbeam trajectories
names = {'ULx','ULy','LLx','LLy',...
         'T1x','T1y','T2x','T2y','T3x','T3y','T4x','T4y',...
         'MNIx','MNIy','MNMx','MNMy'};
for i=I(find(I>0))
  subplot('position',[left_offset+leftPanel_width+sep+0.025 ...
    (nSubplot-k)*height+bottom_offset ...
    rightPanel_width height-sep]);
  plot((1:size(xray(X.m,:),1))'/X.fs,xray(X.m,i),'w'); 
  yAxisLim = [yAxisLim; get(gca,'YLim')];
  plot((1:size(xray,1))'/X.fs,xray(:,i),'b'); hold on;
  %axis([0 max_t yAxisLim(end,:)]);
  ylabel(names(i),'FontSize',13,'Position',[-max_t/20 mean(get(gca,'YLim'))]);
  if k==1
    set(gca,'XAxisLocation','top','Box','on','FontSize',10);
  else
    set(gca,'XTickLabel',[],'Box','on');
  end
  k = k + 1;
end

if ismember(0,I)
  g = subplot('position',[left_offset-0.025 bottom_offset ...
        leftPanel_width 1-top_offset-bottom_offset-0.05]); 
  text(-95,30,['Time: 0 sec'],'FontSize',14); 
  h = axes('position',[left_offset+leftPanel_width+sep+0.025 ...
      bottom_offset rightPanel_width nSubplot*height-sep]);  
  plot([0 0],[-1 1],'r-');
  xlim([0 max_t]);
  axis off
  
  % Clean up memory
  clear X A L 
  
  % Create the QuickTime movie
  if ~isempty(mov)
    MakeQTMovie('start',mov);
    MakeQTMovie('quality',1);
  end
   
  for i=1:mr:size(xray,1)
    % plot vocal tract animation on the left panel
    delete(g);
    g = subplot('position',[left_offset-0.025 bottom_offset ...
      leftPanel_width 1-top_offset-bottom_offset-0.05]); hold on;
    title('Wisconsin X-ray microbeam: Midsagittal Vocal Tract Animation',...
          'FontSize',14);
    text(-95,30,['Time: ' num2str(i/fs_xray,'%1.3f') ' sec'],'FontSize',14);
    XRMBframe(xray(i,:),axis_xy,pal,pha);         
    % plot process bar on the right panel
    delete(h);
    h = axes('position',[left_offset+leftPanel_width+sep+0.025 ...
      bottom_offset rightPanel_width nSubplot*height-sep]);  
    plot([(i-1)/fs_xray (i-1)/fs_xray],[-1 1],'r-');
    xlim([0 max_t]);
    axis off
    % Add frame to the video
    if ~isempty(mov) MakeQTMovie('addfigure'); end
  end
  % Complete the movie by adding the audio
  if ~isempty(mov)
    MakeQTMovie('framerate',ceil(size(xray,1)/mr/max_t));
    MakeQTMovie('addsound',y,fs_wav);
    MakeQTMovie('finish');
  end
  
end

