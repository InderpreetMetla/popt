function grid_analysis = GridHistory(solution,info,iphase,max_error,...
    nodes,tau)
%**********************************************************************%
% This function stores the results and ipopt solution info for each
% grid iteration that the trajectory optimisation problem is solved on.
%
% Inputs:
%
% - solution  : Current solution information in optimal control form
% - info      : Ipopt info (includes number of iterations, cpu secs, etc)
% - iphase    : Phase number
% - max_error : Maximum Relative Error in iphase
% - nodes     : Total number of nodes in iphase
% - tau       : Locations of the nodes in the domain [-1,1]
%
% Outputs:
%
% - grid_analysis : Contains following fields:
%                     - cost    :   Objective for solution
%                     - phase   :   Solution info per phase
%                     - ipopt   :   Info outputted by ipopt
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%**********************************************************************%

grid_analysis.cost = solution.cost;
grid_analysis.phase(iphase).time = solution.phase(iphase).time;
grid_analysis.phase(iphase).state = solution.phase(iphase).state;
grid_analysis.phase(iphase).control = solution.phase(iphase).control;
grid_analysis.phase(iphase).max_relative_error = max_error;
grid_analysis.phase(iphase).number_of_nodes = sum(nodes);
grid_analysis.phase(iphase).node_locations = tau;
grid_analysis.phase(iphase).number_of_intervals = size(nodes,2);

grid_analysis.ipopt.cpu_time = info.cpu;
grid_analysis.ipopt.nlp_iters = info.iter;
grid_analysis.ipopt.eval = info.eval;
grid_analysis.ipopt.status = info.status;

if grid_analysis.ipopt.status == 0
    grid_analysis.ipopt.status_message = ...
        'Optimal Solution Found.';
elseif grid_analysis.ipopt.status == 1
    grid_analysis.ipopt.status_message = ...
        'Solved to an Acceptable Level.';
elseif grid_analysis.ipopt.status == 2
    grid_analysis.ipopt.status_message = ...
        'Converged to a point of local infeasibility. Problem may be infeasible.';
elseif grid_analysis.ipopt.status == 3
    grid_analysis.ipopt.status_message = ...
        'Search Direction Becomes Too Small.';
elseif grid_analysis.ipopt.status == 4
    grid_analysis.ipopt.status_message = ...
        'Diverging Iterates.';
elseif grid_analysis.ipopt.status == 5
    grid_analysis.ipopt.status_message = ...
        'User Requested Stop.';
elseif grid_analysis.ipopt.status == -1
    grid_analysis.ipopt.status_message = ...
        'Maximum Number of Iterations Exceeded.';
elseif grid_analysis.ipopt.status == -2
    grid_analysis.ipopt.status_message = ...
        'Restoration Phase Failed.';
elseif grid_analysis.ipopt.status == -3
    grid_analysis.ipopt.status_message = ...
        'Error in Step Computation.';
elseif grid_analysis.ipopt.status == -10
    grid_analysis.ipopt.status_message = ...
        'Not Enough Degrees of Freedom.';
elseif grid_analysis.ipopt.status == -11
    grid_analysis.ipopt.status_message = ...
        'Invalid Problem Definition.';
elseif grid_analysis.ipopt.status == -12
    grid_analysis.ipopt.status_message = ...
        'Invalid Option.';
elseif grid_analysis.ipopt.status == -13
    grid_analysis.ipopt.status_message = ...
        'Invalid Number Detected.';
elseif grid_analysis.ipopt.status == -100
    grid_analysis.ipopt.status_message = ...
        'Unrecoverable Exception.';
elseif grid_analysis.ipopt.status == -101
    grid_analysis.ipopt.status_message = ...
        'Non-Ipopt Exception Thrown.';
elseif grid_analysis.ipopt.status == -102
    grid_analysis.ipopt.status_message = ...
        'Insufficient Memory.';
elseif grid_analysis.ipopt.status == -199
    grid_analysis.ipopt.status_message = ...
        'Internal Error.';
else
    grid_analysis.ipopt.status_message = ...
        'Unknown';
end

end