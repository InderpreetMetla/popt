function problem = GridRefinement(problem)
%**********************************************************************%
% This function is the grid refinment algorithm.
% The error on the grid is calculated and then it is decided whether
% to increase the number of collocation points in an interval or to 
% sub-divide it. The choice is made in IncreasePolyOrder_or_SubDivide.m
%
% Inputs:
%
% - problem struct
%
% Outputs:
%
% - Updated problem struct with new grid information (nodes, break pts). 
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%**********************************************************************%

solution = problem.solution;
info = problem.info;
iter = problem.grid.iter+1;
auxdata = problem.auxdata;
nphases = problem.nphases;
tol = problem.grid.tol;
End = zeros(1,nphases);
fprintf('\n');
disp('********************************************************************************')
fprintf('Analysis of Grid Iteration %i:',problem.grid.iter+1);
fprintf('\n');
for iphase = 1:nphases
    TauOld = [problem.LGColloc.phase(iphase).Points;1];
    NodesPerInt = problem.grid.phase(iphase).nodes.PerInterval;  
    nInts = size(NodesPerInt,2);
    TemporaryGrid = repmat(struct('NodesPerInt',[],'BreakPts',[]),1,nInts);
    NodesPerIntError = NodesPerInt + 1;
    Idx = [1 cumsum(NodesPerIntError)+ones(1,nInts)];
    BreakPts = problem.grid.phase(iphase).BreakPts;
    Nmin = problem.grid.phase(iphase).nodes.lb;
    Nmax = problem.grid.phase(iphase).nodes.ub;
    % Compute the Relative Errror 
    RelativeError = computeError(solution,auxdata,TauOld,...
        NodesPerInt,BreakPts,nInts,iphase,problem.funcs.Dynamics,...
        problem.scaling.stateScales);
    MRE = zeros(1,nInts);
    for iInt = 1:nInts
        MRE(iInt) = MaxRelError(RelativeError,Idx,iInt);
    end
    fprintf('    Max Relative Error in Phase %i is %i',iphase,max(MRE));
    fprintf('\n');

    ErrorSatisfactionFlag = (RelativeError <= tol);
    ErrorSatisfactionFlag = all(ErrorSatisfactionFlag(:));

    if ErrorSatisfactionFlag == 0
        for iInt = 1:nInts
            TemporaryGrid = ...
                IncreasePolyOrder_or_SubDivide(NodesPerInt,BreakPts,...
                TemporaryGrid,MRE,iphase,iInt,Nmin,Nmax,tol);
        end
        NodesPerIntNew = [];
        BreakPtsNew = [];
        if ~isempty(TemporaryGrid)
            for iInt = 1:nInts
                if ~isempty(TemporaryGrid(iInt).BreakPts)
                    BreakPtsNew = [BreakPtsNew,...
                        TemporaryGrid(iInt).BreakPts];
                else
                    BreakPtsNew = [BreakPtsNew,...
                        BreakPts(iInt:iInt+1)];
                end
                if ~isempty(TemporaryGrid(iInt).NodesPerInt)
                    NodesPerIntNew = [NodesPerIntNew,...
                        TemporaryGrid(iInt).NodesPerInt];
                else
                    NodesPerIntNew = [NodesPerIntNew,...
                        NodesPerInt(iInt)];
                end
            end
        else
            NodesPerIntNew = NodesPerInt;
            BreakPtsNew    = BreakPts;
        end
        ExitFlag = 0;
    else
        NodesPerIntNew = NodesPerInt;
        BreakPtsNew    = BreakPts;
        ExitFlag = 0;
        fprintf('    Error Tolerance Satisfied in Phase %i.',iphase);
        fprintf('\n');
    end
    
    if isequal(NodesPerIntNew,NodesPerInt) && ...
            isequal(BreakPtsNew,BreakPts) && ...
            isequal(ExitFlag,0) && (isequal(info.status,0) || ...
            isequal(info.status,1))
        End(iphase) = 0;
    else
        End(iphase) = 1;
    end
    
    problem.grid.phase(iphase).nodes.PerInterval = NodesPerIntNew;
    problem.grid.phase(iphase).BreakPts = unique(BreakPtsNew);
    
    % Store an analysis of the Grid along with collocation history
    Grid_Analysis = GridHistory(solution,info,iphase,...
        max(MRE),NodesPerInt,TauOld);
end
problem.grid_analysis.iteration(iter) = Grid_Analysis;
if any(End == 1)
    problem.ErrorFlag = 1;
else
    problem.ErrorFlag = 0;
end
disp('********************************************************************************')
fprintf('\n');
end