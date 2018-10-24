function [solution, plotdata] = popt(problem)
%*******************************************************************%
%                      POPT: Pseudospectral OPTimiser
%*******************************************************************%
% Distributed under the MIT License
% License information is available in LICENSE.txt
%
% Copyright (c) 2018 Inderpreet Metla
%*******************************************************************%
%
% Uses Radau Pseudospectral method to solve trajectory optimisation
% problems. 
%
%*******************************************************************%
disp('*************************************************************')
disp('*                                                           *')
disp('*               POPT: Pseudospectral OPTimiser              *')
disp('*                 Written by Inderpreet Metla               *')
disp('*                                                           *')
disp('*                         Version 1.1                       *')
disp('*                                                           *')
disp('*************************************************************')
disp(' ')
%*******************************************************************%
% Assign Binary Mex file-name extension for current OS
%*******************************************************************%
ext = mexext;

%*******************************************************************%
% Check for Ipopt installation
%*******************************************************************%
if isempty(which(strcat('ipopt.',ext)))
  error('Ipopt not found or is not installed correctly')
end

%*******************************************************************%
% Input Validation
%*******************************************************************%
inputValidator(problem);

%*******************************************************************%
% Initial Grid Setup
%*******************************************************************%
problem = GridSetup(problem); 

%*******************************************************************%
% Determine the Number of Decision Variables and Constraints 
%*******************************************************************%
problem = NumofDecVarsandCsts(problem);

%*******************************************************************%
% Create a dummy auxdata variable to avoid errors
%*******************************************************************%
problem.auxdata.dummy = 0;

%*******************************************************************%
% Dependencies
%*******************************************************************%
problem = ProblemDependencies(problem);

%*******************************************************************%
% Solver and Grid Refinement Loop
%*******************************************************************%
ErrorFlag = 1;
problem.cpu_time.grids = [];
problem.cpu_time.total = 0;
while ErrorFlag == 1 && ...
        (problem.grid.iter <= problem.grid.max_refine)
%*******************************************************************%
% Determine the Number of Nodes in each phase of the problem
%*******************************************************************%
problem = NumofNodes(problem);

%*******************************************************************%
% Extract Guess
%*******************************************************************%
problem = InitialGuessExtractor(problem);

%*******************************************************************%
% Determine the bounds on the decision variables and constraints
%*******************************************************************%
problem = DecVarAndCstBounds(problem);

%*******************************************************************%
% Scale the NLP
%*******************************************************************%
problem = ScaleNlp(problem);

%*******************************************************************%
% IPOPT Options -Defaults and User Defined
%*******************************************************************%
problem = optionSetup(problem);

%*******************************************************************%
% Initial Guess (scaled) 
%*******************************************************************%
x0 = problem.DecVar0;
x0 = x0.*problem.scaling.decvar_scales + ...
    problem.scaling.decvar_shifts;

%*******************************************************************%
% Assign problem to a Global Variable called MAIN
%*******************************************************************%
global MAIN
MAIN = problem;

%*******************************************************************%
% Validate that Objective Function and Constraint Function are 
% correctly set up
%*******************************************************************%
objFunValidator(x0,problem.auxdata);
cstFunValidator(x0,problem.auxdata);

%*******************************************************************%
% Derivative Method
%*******************************************************************%
funcs = DerivMethod(problem);

%*******************************************************************%
% Display Grid Iteration Number
%*******************************************************************%
disp('-------------------------------------------------------------')
fprintf('Solving on Grid Iteration %i ',problem.grid.iter+1)
fprintf('\n')
for iphase = 1:problem.nphases
fprintf('	Phase %i:',iphase)   
fprintf('\n')
    fprintf('       Total Number of Collocation Points = %i', ...
        problem.nodes(iphase))
    fprintf('\n')
    fprintf('       Total Number of Grid Intervals     = %i',...
        length(problem.grid.phase(iphase).nodes.PerInterval))
    fprintf('\n')
end
disp('-------------------------------------------------------------')

%*******************************************************************%
% Call IPOPT
%*******************************************************************%
[x, info] = ipopt(x0,funcs,problem.options);

%*******************************************************************%
% Solution Post Processing
% Convert the solution for the Non Linear Program back to 
% Trajectory Optimisation form
%*******************************************************************%
problem.result = (x-problem.scaling.decvar_shifts)./...
    problem.scaling.decvar_scales; 
problem.info   = info;
problem.cpu_time.grids = [problem.cpu_time.grids,info.cpu];
problem.cpu_time.total = problem.cpu_time.total+info.cpu;
if nargout == 2
    [solution, plotdata] = SolutionUnpacker(problem);
    PrintFlag = 1;
else
    solution = SolutionUnpacker(problem);
    PrintFlag = 2;
end

problem.solution = solution;
problem = GridRefinement(problem);
ErrorFlag = problem.ErrorFlag;
problem.grid.iter = problem.grid.iter + 1;
problem.guess = problem.solution;
clear problem.result
% End grid refinement while loop:
end 

% Append Grid Analysis History to the solution
solution.grid_analysis = problem.grid_analysis;

if problem.grid.iter > problem.grid.max_refine
    disp('Max Number of Grid Refinements Exhausted.')
    disp('Problem May Not Be Solved.')
end

outputdisplay(PrintFlag)
end