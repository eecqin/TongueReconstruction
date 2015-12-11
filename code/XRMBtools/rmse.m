% [E,Estd] = rmse(A,B)
% Root-mean-square error between rowwise data sets A and B
%
% In:
%   A: NxD matrix of estimated articulatory trajectories.
%   B: NxD matrix of true articulatory trajectories.
% Out:
%   E: Dx1 vector of RMS error for each articulatory channel.
%   Estd: Dx1 vector of RMS error standard deviation for each articulatory
%      channel.
% 
% Copyright (c) 2008 by Chao Qin and Miguel A. Carreira-Perpinan

function [E,Estd] = rmse(A,B)
  
% Vectorised version to compute root-mean-square error
C = A-B; E = sqrt(mean(C.^2,1))'; Estd = sqrt(std(C.^2,0,1))';

