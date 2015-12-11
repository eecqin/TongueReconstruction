% [E,A,L,L2] = MOCHAread(uttfile) 
% 
% Read EMA data, acoustic wave and phonetic labels associated with an
% utterance. Expects the corresponding 3 files (.ema, .wav, .lab) to exist.
% 
% In: 
%  uttfile: utterance filename.
% Out:
%  E,A,L: see MOCHAplot.
%  L2: phonetic labels for each frame.

% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2009 by Chao Qin and Miguel A. Carreira-Perpinan

function [E,A,L,L2] = MOCHAread(uttfile)

% Load EMA data in ASCII format
if exist([uttfile '.ema'],'file')
  E.dat = load([uttfile '.ema'],'-ascii');
  % Re-order and normalize orginal EMA data
  E.dat = E.dat(:,[4 9 5 10 6 11 7 12 13 18 14 19 15 20 3 8 17 22])/100;
  E.fs = 500; 
else
  E = []; fprintf(1,[uttfile '.ema does not exist!']);
end

% Load acoustic wave
if exist([uttfile '.wav'],'file')
  [y,fs,nbits] = wavread([uttfile '.wav']);
  A.dat = y; A.fs = fs; A.nbits = nbits;
  if size(E.dat,1)/E.fs<=size(A.dat,1)/A.fs
    max_t = size(E.dat,1)/E.fs; 
  else
    max_t = size(A.dat,1)/A.fs;
  end
else
  A = []; fprintf(1,[uttfile '.wav does not exist!']);
end

% Load phonetic label
if exist([uttfile '.lab'],'file')
  k = 1; 
  flab = fopen([uttfile '.lab'],'r');
  while ~feof(flab)
    tline = fgetl(flab);
    [ts,tmp] = strtok(tline); [te,tmp] = strtok(tmp);
    [phone,tmp] = strtok(tmp);
    if strcmp(phone,'breath') phone = 'br'; end;
    te = num2str(min(str2num(te),max_t));  
    L(k).phone = phone; L(k).startt = ts; L(k).endt = te;         
    for i=round(str2num(ts)*100)+1:round(str2num(te)*100) 
      L2{i} = phone; 
    end;
    k = k + 1;
  end
  fclose(flab);
else
  L = []; L2 = []; fprintf(1,[uttfile '.lab does not exist!']);
end

