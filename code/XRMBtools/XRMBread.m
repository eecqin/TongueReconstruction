% [X,A,L,L2] = XRMBread(uttfile) 
% 
% Read X-ray microbeam data, acoustic wave and phonetic labels associated with
% an utterance. Expects the corresponding 3 files (.ema, .wav, .lab) to exist.
% 
% In: 
%  uttfile: utterance filename.
% Out:
%  X,A,L: see XRMBplot.
%  L2: phonetic labels for each frame.

% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2008 by Chao Qin and Miguel A. Carreira-Perpinan

function [X,A,L,L2] = XRMBread(uttfile)

% Load X-ray microbeam data in ASCII format
if exist([uttfile '.txy'],'file')
  X.dat = load([uttfile '.txy'],'-ascii');
  % select only XRAY data
  X.dat = X.dat(:,2:end)/1000;   
  X.fs = 145.6; 
  X.m = [];
  for i=1:16
    X.m = [X.m; find(X.dat(:,i)==1000)];
    X.dat(find(X.dat(:,i)==1000),i) = NaN;
  end
  X.m = unique(X.m); X.m = sort(X.m);
  X.m = setdiff([1:size(X.dat,1)]',X.m);
else
  X = [];
  fprintf(1,[uttfile '.txy does not exist!\n']);
end

% Load acoustic wave
if exist([uttfile '.wav'],'file')
  [y,fs,nbits] = wavread([uttfile '.wav']);
  A.dat = y; A.fs = fs; A.nbits = nbits;
  if size(X.dat,1)/X.fs<=size(A.dat,1)/A.fs
    max_t = size(X.dat,1)/X.fs; 
  else
    max_t = size(A.dat,1)/A.fs;
  end
else
  A = [];
  fprintf(1,[uttfile '.wav does not exist!\n']);
end

% Load phonetic label
if exist([uttfile '.lab'],'file')
  k = 1; 
  flab = fopen([uttfile '.lab'],'r');
  while ~feof(flab)
    tline = fgetl(flab);
    [ts,tmp] = strtok(tline);
    [te,tmp] = strtok(tmp);
    [phone,tmp] = strtok(tmp);
    if strcmp(phone,'breath') phone = 'br'; end;
    te = num2str(min(str2num(te),max_t));  
    L(k).phone = phone; L(k).startt = ts; L(k).endt = te;
    for i=round(str2num(ts)*X.fs)+1:round(str2num(te)*X.fs) 
      L2{i} = phone; 
    end;
    k = k + 1;
  end
  fclose(flab);
else
  L = []; L2 = []; fprintf(1,[uttfile '.lab does not exist!\n']);
end

