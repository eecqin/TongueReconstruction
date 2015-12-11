% Driver file to test functions in MOCHAtools

% Auxiliary functions needed by XRMBplot from our speech_analysis toolbox
% (assumed installed in ../speech_analysis)

% Indices of articulators to display
I = [1:8];		% Display all articulators
axis_xy = [-20 80 -40 20];

% Test MOCHAread
[E,A,L] = MOCHAread('dat/fsew0_001');

% Test MOCHAscat (its display looks similar to MOCHAtrace's because the
% input EMA data is from only one utterance).
figure(1); clf; MOCHAscat(E.dat,I,axis_xy);

% Test MOCHAtrace
figure(2); clf; MOCHAtrace(E.dat,I,axis_xy);

% Test MOCHAframe
for i=1:size(E.dat,1)
  figure(3); clf; hold on;
  text(35,60,['Time: ' num2str(i/E.fs,'%1.3f') ' sec'],'FontSize',14);
  MOCHAframe(E.dat(i,:),I,axis_xy);
end

% Test MOCHAdownsamp
wn = 0.8; order = 100; downsamp_type = 'filtfilt';
E2.fs = 25;		% Downsample EMA from 500Hz to 25Hz
E2.dat = MOCHAdownsamp(E.dat,round(E.fs/E2.fs),downsamp_type,wn,order);

% Test MOCHAplot and create a movie 'fsew0_001.mov'.
I = [0 -1 -2 -3 2:2:14];	% Display selected acoustics and articulators
mr = 40;			% Moving rate for the process bar
figure(4); clf; MOCHAplot(E,A,L,I,[],[],[],[],['fsew0_001.mov'],mr);

