function problem  = GridSetup(problem)
%******************************************************************%
% Function sets up the fields for the initial grid that the problem
% will be solved on. Sets default values if not user defined.
% This includes the grid error tolerance, maximum number of grid
% refinements, upper and lower bounds on the number of nodes per
% grid interval as the number of nodes per interval in the initial
% grid. If the nodes in the initial grid are supplied as a scalar 
% value, then there is a single interval. If a vector of nodes is 
% provided, then the initial grid is broken into multiple intervals.
% 
% Inputs:
%
% - problem struct
%
% Outputs:
% - problem struct with:
%       - .grid.tol        : Error Tolerance for the grid in 
%                            every phase
%       - .grid.max_refine : Maximum number of grid refinements
%       - .grid.iter       : Grid iteration
%       - .grid.phase:
%           - nodes.initialgrid : Node distribution in initial grid
%           - nodes.lb          : Lower Bound on nodes for the 
%                                 grid refinement algorithm
%           - nodes.ub          : Upper Bound on nodes for the 
%                                 grid refinement algorithm
%           - nodes.BreakPts    : Where the intervals will be broken
%                                 along the [-1,1] Radau time domain
%           - nodes.PerInterval : Same as nodes.initialgrid expect 
%                                 this will update through the
%                                 solution as the grid is refined 
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%******************************************************************%

% Number of phases
nphases = length(problem.bounds.phase);

if isfield(problem,'grid')
    if ~isfield(problem.grid,'tol') ||...
            isempty(problem.grid.tol)
        problem.grid.tol = 1e-7;
    end
    if ~isfield(problem.grid,'max_refine') ||...
            isempty(problem.grid.max_refine)
        problem.grid.max_refine = 10;
    end
    if isfield(problem.grid,'phase')
        Grid = problem.grid.phase;
        for iphase = 1:nphases
            % Grid(iphase).BreakPts = [-1,1];
            if ~isfield(Grid(iphase),'nodes') || ...
                    isempty(Grid(iphase).nodes)
                Grid(iphase).nodes.initialgrid = 4*ones(1,10);
                Grid(iphase).nodes.lb = 3;
                Grid(iphase).nodes.ub = 10;
            else
                if ~isfield(Grid(iphase).nodes,'initialgrid') ||...
                        isempty(Grid(iphase).nodes.initialgrid)
                    Grid(iphase).nodes.initialgrid = 4*ones(1,10);
                end
                if ~isfield(Grid(iphase).nodes,'lb') ||...
                        isempty(Grid(iphase).nodes.lb)
                    Grid(iphase).nodes.lb = 3;
                end
                if ~isfield(Grid(iphase).nodes,'ub') ||...
                        isempty(Grid(iphase).nodes.ub)
                    Grid(iphase).nodes.ub = 10;
                end
            end
            nInts = length(Grid(iphase).nodes.initialgrid);
            Grid(iphase).BreakPts = -1:2/nInts:1;
        end
        problem.grid.phase = Grid;
    else
        for iphase = 1:nphases
            problem.grid.phase(iphase).nodes.initialgrid = 4*ones(1,10);
            problem.grid.phase(iphase).nodes.lb = 3;
            problem.grid.phase(iphase).nodes.ub = 10;
            problem.grid.phase(iphase).BreakPts = -1:2/10:1;
        end
    end
else
    problem.grid.tol = 1e-7;
    problem.grid.max_refine = 10;
    for iphase = 1:nphases
        problem.grid.phase(iphase).nodes.initialgrid = 4*ones(1,10);
        problem.grid.phase(iphase).nodes.lb = 3;
        problem.grid.phase(iphase).nodes.ub = 10;
        problem.grid.phase(iphase).BreakPts = -1:2/10:1;
    end
end
problem.grid.iter = 0;

for iphase = 1:nphases
    problem.grid.phase(iphase).nodes.PerInterval = ...
        problem.grid.phase(iphase).nodes.initialgrid;
end

end