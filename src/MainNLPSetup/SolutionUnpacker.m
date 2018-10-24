function [solution, plotdata] = SolutionUnpacker(problem)
%******************************************************************%
% This function computes two things. 
%
% (1) It takes the single vector  solution of the optimisied 
% decision variables and converts it back into the form of  
% a trajectory optimisation problem. i.e it converts the
% non-linear program solution to optimal control form.
%
% (2) It uses barycentric interpolation to construct an output 
% struct called plotdata that can be used to quickly form plots. 
% plotdata is only created if requested by the user. 
% 
% Inputs:
% - problem :	User supplied problem struct + new fields from 
%               previous functions plus the result from ipopt
%
% Outputs:
%  - solution:	struct with following fields
%
%       - solution.phase        :   struct with fields of .time,
%                                   .state, .control which are 
%                                   the solutions at the GPM
%                                   collocation points
%       - solution.cost         :	Total objective cost
%       - solution.bndcost      : 	Boundary (Mayer) objective cost
%       - solution.pathcost     :	Path (Lagrangian) objective cost
%       - solution.status       :   ipopt exit condition
%       - solution.num_iters    :   Number of ipopt iterations
%       - solution.f_evals      :   Number of ipopt function evals
%       - solution.cpu_time 	:   Ipopt cpu time
%
%  - plotdata:	struct with following fields
%
%       - plotdata.phase.time    :	Time data struct
%       - plotdata.phase.state   : 	State data struct
%       - plotdata.phase.control :	Control data struct
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%******************************************************************%

% Call NLP Solution
x           = problem.result;
% Struct of Indexes of the State (x), Control (u) and Time (t) 
Idx_xut     = problem.Idx;
% Number of Nodes 
nodes       = problem.nodes;
% Number of Phases
nphases     = problem.nphases;
% Matrix of number of States, Controls, Path and Bnd Csts
sizes       = problem.sizes;
% Number of States in Each Phase (row vector)
nStates     = sizes(1,:);
% Number of Controls in Each Phase (row vector)
nControls   = sizes(2,:);
% GPM Collocation Information
LGColloc    = problem.LGColloc;
%******************************************************************%
BndCost  = 0;
PathCost = 0;
Cost     = 0;
auxdata = problem.auxdata;
for iphase = 1:nphases
    states      = x(Idx_xut.phase(iphase).stateIdx);
    controls    = x(Idx_xut.phase(iphase).controlIdx);
    t0          = x(Idx_xut.phase(iphase).timeIdx(1));
    tf          = x(Idx_xut.phase(iphase).timeIdx(2));
    LGPoints    = LGColloc.phase(iphase).Points;
    Tau         = [LGPoints;1];
    Time        = 0.5 * (tf - t0) .* (Tau + 1) + t0;
    TimeGR      = Time(1:end-1);
    
    solution.phase(iphase).time = Time;
    
    solution.phase(iphase).state = reshape(states,nodes(iphase)+1,...
        nStates(iphase));

    if ~isempty(controls)
        controlGRMatrix = reshape(controls,nodes(iphase),...
            nControls(iphase));
        solution.phase(iphase).control  = interp1(TimeGR,...
            controlGRMatrix,Time,'pchip','extrap');
    else
        controlGRMatrix = [];
        solution.phase(iphase).control  = [];
    end
    
    if nargout == 2
        %********************************************************%
        % Plotting Solution
        %********************************************************%
        if nodes(iphase) < 30
            N = 100;
        elseif nodes(iphase) >= 30 && nodes(iphase) < 80
            N = 120;
        elseif nodes(iphase) >= 80 && nodes(iphase) < 150
            N = 200;
        elseif nodes(iphase) >= 150
            N = min(nodes(iphase),200);
        end
        plotdata.phase(iphase).time    = linspace(t0,tf,N).';
        plotdata.phase(iphase).state   = ...
            interp1(Time,solution.phase(iphase).state,...
            plotdata.phase(iphase).time,'pchip');
%             barylag([time,solution.phase(iphase).state],...
%             plotdata.phase(iphase).time);
        if ~isempty(controls)
            plotdata.phase(iphase).control = ...
                interp1(Time,solution.phase(iphase).control,...
                plotdata.phase(iphase).time,'pchip');
%                 barylag([time,solution.phase(iphase).control],...
%                 plotdata.phase(iphase).time);
        else
            plotdata.phase(iphase).control = [];
        end
    end
    %********************************************************%
    % Cost Function Evaluation Intialisation
    %********************************************************%
    x0 = solution.phase(iphase).state(1,:);
    xf = solution.phase(iphase).state(end,:);
    stateGRMatrix = solution.phase(iphase).state(1:end-1,:);

%********************************************************%
% Cost Calculation
%********************************************************%
    auxdata.iphase = iphase;
    Lagrangian = problem.funcs.PathObj(TimeGR,stateGRMatrix,...
        controlGRMatrix,auxdata);
    Mayer  = problem.funcs.BndObj(t0,tf,x0,xf,auxdata);
	  
    LGWeights   = LGColloc.phase(iphase).Weights;
    Quadrature  = 0.5 * (tf - t0) * LGWeights * Lagrangian;
    Cost        = Cost + Mayer + Quadrature;
    BndCost     = BndCost + Mayer;
    PathCost    = PathCost + Quadrature;
end

solution.cost       = Cost;
solution.bndcost	= BndCost;
solution.pathcost	= PathCost;

solution.status         = problem.info.status;
solution.nlp_iters      = problem.info.iter;
solution.f_evals        = problem.info.eval;
solution.cpu_time.grids = problem.cpu_time.grids;
solution.cpu_time.total = problem.cpu_time.total;
solution.grid_iters     = problem.grid.iter+1;
end