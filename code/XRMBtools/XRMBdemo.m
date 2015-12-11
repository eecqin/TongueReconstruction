% Driver file to test functions in XRMBtools

% Auxiliary functions needed by XRMBplot from our speech_analysis toolbox
% (assumed installed in ../speech_analysis)
addpath /home/cqin/matlab_files/src/speech/cqin

% Test XRMBread
[X,A,L] = XRMBread('dat/jw11_tp105');
load('dat/jw11_pha.dat','-ascii'); % pharyngeal wall
load('dat/jw11_pal.dat','-ascii'); % palate outline 

% Test XRMBscat (its display looks similar to XRMBtrace's because the
% input X-ray microbeam data is from only one utterance).
figure(1); clf; XRMBscat(X.dat(X.m,:),[],jw11_pal/1000,jw11_pha/1000); 

% Test XRMBtrace
I = [1:8];		% Display all articulators
figure(2); clf; XRMBtrace(X.dat(X.m,:),I,[],jw11_pal/1000,jw11_pha/1000);

% Test XRMBframe
for i=1:size(X.dat,1)
  figure(3); clf; hold on;
  text(-95,30,['Time: ' num2str(i/X.fs,'%1.3f') ' sec'],'FontSize',14); 
  XRMBframe(X.dat(i,:),[],jw11_pal/1000,jw11_pha/1000);
end

% Test XRMBdownsamp
% Important note: XRMBdownsamp works only with "good" X-ray microbeam data.
% In the XRMB database, bad values are coded as 1000, and they often occur
% and the beginning or end of each utterance.
wn = 0.8; order = 100; downsamp_type = 'filtfilt';
X2.fs = 100;		% Downsample X-ray microbeam from 146Hz to 100Hz
X2.dat = XRMBdownsamp(X.dat,round(X.fs/X2.fs),downsamp_type,wn,order);
X2.m = [1:size(X2.dat,1)]';

% Test XRMBplot and create a movie 'jw11_tp105.mov'.
% You will see occasionally missing articulators in the midsagittal plane
% during the animation, caused by bad data (code 1000 in the XRMB database).
I = [0 -1 -2 -3 -4 1:2:13];	% Display selected acoustics and articulators
mr = 32;			% Moving rate for the process bar
figure(4); clf;
XRMBplot(X,A,L,I,[],jw11_pal/1000,jw11_pha/1000,[],['jw11_tp105.mov'],mr);


spk = 'jw11';
fid = fopen(['/home/cqin/db/XRMB/processed_ubdb/' spk '/' spk '_2.list'],'r');
figure(5);
while feof(fid)~=1
  filename = fgetl(fid)
  [X,A,L] = XRMBread(['/home/cqin/db/XRMB/processed_ubdb/' spk '/' filename]); 
  
  set(0,'CurrentFigure',5); clf; XRMBplot(X,A,L,[1:1:16]);drawnow;pause
end
fclose(fid);
