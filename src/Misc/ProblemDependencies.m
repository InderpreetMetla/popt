function problem = ProblemDependencies(problem)
%**********************************************************************%
% Function computes the variable dependencies for the optimal control 
% functions. It uses the NaN approach. If a variable is replaced with 
% NaN and a function is dependent on that variable, matlab will output
% NaN. This allows one to determine which variables a multivariable
% function is dependent on. This code only runs for states, controls,
% and time as single units and does not check for dependency of 
% individual nodes. This means, if f(z) is a dynamics function and the 
% problem has 3 states, if f(z) is dependent on state 1, 2 or 3 it is 
% assumed to be dependent on all 3 at all node locations. 
% Similarly for controls. 
%
% Inputs:
%
% - problem struct
%
% Outputs:
% - problem struct with:
%       - .dependency  : Contains the dependency of the dynamics and 
%                        path constraints function on time, 
%                        state and control. Also containts the 
%                        dependency of the boundary constraint function
%                        on initial/final time and state.
%                        This is outputted for every phase. 
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%**********************************************************************%

% Number of Phases
nphases = problem.nphases;
% Matrix of number of States, Controls, Path and Bnd Csts
sizes = problem.sizes;
% Number of States in Each Phase (row vector)
nStates = sizes(1,:);
% Number of Controls in Each Phase (row vector)
nControls = sizes(2,:);
for iphase = 1:nphases
    Nx = nStates(iphase);
    Nu = nControls(iphase);
    auxdata = problem.auxdata;
    auxdata.iphase = iphase;
    %*************************************************************%
    % Check if Dynamics and Path Constraints Dependent on State
    %*************************************************************%
    state = NaN*ones(1,Nx);
    control = rand(1,Nu) + eps;
    time = rand + eps;
    % Dynamics wrt state
    Dyn = problem.funcs.Dynamics(time,state,control,auxdata);
    problem.dependency(iphase).dynamics.state = any(isnan(Dyn));
    % Path Cst wrt state
    PathCst = problem.funcs.PathCst(time,state,control,auxdata);
    problem.dependency(iphase).pathcst.state = any(isnan(PathCst));
    %*************************************************************%
    % Check if Dynamics and Path Constraints Dependent on Control
    %*************************************************************%
    state = rand(1,Nx) + eps;
    control = NaN*ones(1,Nu);
    time = rand + eps;
    % Dynamics wrt control
    Dyn = problem.funcs.Dynamics(time,state,control,auxdata);
    problem.dependency(iphase).dynamics.control = any(isnan(Dyn));
    % Path Cst wrt control
    PathCst = problem.funcs.PathCst(time,state,control,auxdata);
    problem.dependency(iphase).pathcst.control = any(isnan(PathCst));
    %*************************************************************%
    % Check if Dynamics and Path Constraints Dependent on Time
    %*************************************************************%
    state = rand(1,Nx)+ eps;
    control = rand(1,Nu)+ eps;
    time = NaN;
    % Dynamics wrt control
    Dyn = problem.funcs.Dynamics(time,state,control,auxdata);
    problem.dependency(iphase).dynamics.time = any(isnan(Dyn));
    % Path Cst wrt control
    PathCst = problem.funcs.PathCst(time,state,control,auxdata);
    problem.dependency(iphase).pathcst.time = any(isnan(PathCst));
    %*************************************************************%
    % Check if Bnd Constraint Dependent on Inputs
    %*************************************************************%
    % Initial Time
    t0 = NaN;
    tf = rand+ eps;
    x0 = rand(1,Nx) + eps;
    xf = rand(1,Nx) + eps;
    BndCst = problem.funcs.BndCst(t0,tf,x0,xf,auxdata);
    problem.dependency(iphase).bndcst.initialtime = any(isnan(BndCst));
    % Final Time
    t0 = rand+ eps;
    tf = NaN;
    x0 = rand(1,Nx) + eps;
    xf = rand(1,Nx) + eps;
    BndCst = problem.funcs.BndCst(t0,tf,x0,xf,auxdata);
    problem.dependency(iphase).bndcst.finaltime = any(isnan(BndCst));
    % Initial State
    t0 = rand  + eps;
    tf = rand + eps;
    x0 = NaN*ones(1,Nx);
    xf = rand(1,Nx) + eps;
    BndCst = problem.funcs.BndCst(t0,tf,x0,xf,auxdata);
    problem.dependency(iphase).bndcst.initialstate = any(isnan(BndCst));
    % Final State
    t0 = rand + eps;
    tf = rand + eps;
    x0 = rand(1,Nx) + eps;
    xf = NaN*ones(1,Nx);
    BndCst = problem.funcs.BndCst(t0,tf,x0,xf,auxdata);
    problem.dependency(iphase).bndcst.finalstate = any(isnan(BndCst));
end