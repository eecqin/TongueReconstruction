% x2 = XRMBdownsamp(x[,ratio,type,wn,order]) 
%
% Performs downsampling on original X-ray microbeam trajectories from an
% utterance. It allows two types of downsampling, based on the Matlab Signal
% Processing Toolbox functions 'resample' and 'filtfilt'. 'resample' can
% induce artifacts on the boundaries of temporal trajectories but these
% trajectories are more smooth; 'filtfilt' produces natural boundaries but
% less smooth temporal trajectories. Overall, their differences are minor.
% 
% In: 
%  x: see XRMBframe.
%  ratio: ratio of downsampling. Default: 1.45.
%  type: type of downsampling. Default: 'resample'.
%  wn: normalized cut-off frequency for double filtering. Default: 0.8.
%  order: order of the double filter. Default: 100.
% 
% Out:
%  x2: ?x16 vector of downsampled X-ray microbeam trajectories.

% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2008 by Chao Qin and Miguel A. Carreira-Perpinan

function x2 = XRMBdownsamp(x,ratio,type,wn,order)

% --------------- Argument defaults ---------------- %
if ~exist('ratio','var') | isempty(ratio) ratio = 1.45; end;
if ~exist('type','var') | isempty(type) type = 'resample'; end;
if ~exist('wn','var') | isempty(wn) wn = 0.8; end;
if ~exist('order','var') | isempty(order) order = 100; end;
% --------------- Argument defaults ---------------- %

[N,D] = size(x); 
x2 = zeros(ceil(N/ratio),D);

% Downsample X-ray trajectories
switch type
 case 'resample'
  x2 = resample(x,1,ratio);
 case 'filtfilt'
  b = fir1(order,wn); a = 1;
  for i=1:D
    x2(:,i) = downsample(filtfilt(b,a,x(:,i)),ratio);
  end
 otherwise
  error('Wrong downsampling type!');
end

