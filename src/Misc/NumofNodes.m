function problem = NumofNodes(problem)
%********************************************************************%
% Determine the total number of nodes for each phase             
%
% Inputs:
%
% - problem struct
%
% Outputs:
%
% - problem struct with number of nodes in each phase ouputed as row
%   a vector
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%********************************************************************%

% Number of phases
nphases = length(problem.bounds.phase);
% Rename initial grid struct for better readability
Grid = problem.grid.phase;
% Set up problem.nodes field 
problem.nodes = zeros(1,nphases);

% Append Total Number of Nodes to problem.nodes
for iphase=1:nphases
    problem.nodes(iphase) = sum(Grid(iphase).nodes.PerInterval);
end

end