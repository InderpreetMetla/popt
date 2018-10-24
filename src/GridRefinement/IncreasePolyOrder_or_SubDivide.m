function TemporaryGrid = ...
    IncreasePolyOrder_or_SubDivide(NodesPerInt,BreakPts,...
    TemporaryGrid,MRE,iphase,iInt,Nmin,Nmax,tol)
%**********************************************************************%
% This function decides whether to increase the number of collocation 
% points in an interval or to sub-divide it. 
%
% Outputs:
%
% - TemporaryGrid : Contains the information for the new grid. Stored 
%                   in a temporary grid so the original grid struct
%                   doesn't get changed until required.
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%**********************************************************************%
N = NodesPerInt(iInt); % -1
IntStart = BreakPts(iInt);
IntEnd   = BreakPts(iInt+1);
MaxRelError = MRE(iInt);
if MaxRelError > tol
    P = ceil(log(MaxRelError/tol) / log(N));
    if (N+P) <= Nmax
        % Increase the order of the polynomial in this interval and
        % keep the number of segments the same
        TemporaryGrid(iInt).NodesPerInt = N + P;
        TemporaryGrid(iInt).BreakPts    = [IntStart,IntEnd];
        fprintf(strcat('    Increasing Number of Collocation',...
            ' Points in Interval %i from %i to %i'),iInt,N,N+P);
        fprintf('\n');
    else
        % Subdivide this interval and set #Nodes = Nmin in each Division
        % Limit the number of subintervals to 12
        nSubIntervals = min(max(ceil((N+P)/Nmin),2),12);
        GridSplitPoints = IntStart:(IntEnd-IntStart)/nSubIntervals:IntEnd;
        TemporaryGrid(iInt).NodesPerInt = Nmin*ones(1,nSubIntervals);
        TemporaryGrid(iInt).BreakPts    = GridSplitPoints;
        fprintf(strcat('    Refining Interval %i Into %i',...
            ' Sub Intervals'),iInt,nSubIntervals);
        fprintf('\n');
    end
else
    % Keep this interval unchanged
    fprintf(strcat('    Error Tolerance Satisfied',...
            ' in Interval %i of Phase %i'),iInt,iphase);
    fprintf('\n');
end
end