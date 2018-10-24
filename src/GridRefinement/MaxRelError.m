function MRE = MaxRelError(RelativeError,Idx,iInt)
%**********************************************************************%
% This function computes the maximum relative error across every time 
% point of every state
% 
% Inputs: 
%
%   - RelativeError : The Relative Error from computeError.m
%   - Idx           : The index for the relative errors in the 
%                     current interval
%   - iInt          : Current interval 
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%**********************************************************************%

IntervalRelError = RelativeError(Idx(iInt):Idx(iInt+1),:);
MRE = max(max(IntervalRelError));
end