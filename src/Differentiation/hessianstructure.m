function Htotal = hessianstructure()
%*************************************************************%
% This function computes the sparsity structure of the Hessian 
%
% Outputs: 
%    - Htotal  : Lagrangian Hessian Sparsity Pattern
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%*************************************************************%
global MAIN

% Number of Nodes 
Nodes = MAIN.nodes;
% Matrix of number of States, Controls, Path and Bnd Csts
sizes = MAIN.sizes;
% Number of States in Each Phase (row vector)
nStates = sizes(1,:);
% Number of Controls in Each Phase (row vector)
nControls = sizes(2,:);
% Allocate Memory for Total Jacobian
Htotal = sparse(size(MAIN.DecVar0,1),size(MAIN.DecVar0,1));
shift = 0;
for iphase = 1:MAIN.nphases
    N   = Nodes(iphase);
    Nx  = nStates(iphase);
    Nu  = nControls(iphase);
    
    StateStateBlock = sparse(diag(ones(N+1,1)));
    StateStateBlock(1,end)=1;
    StateStateBlock(end,1)=1; 
    StateStateBlock     = repmat(StateStateBlock,Nx,Nx);
    ControlStateBlock   = repmat(sparse([diag(ones(N,1));zeros(1,N)]),Nx,Nu);
    StateControlBlock   = ControlStateBlock.';
    ControlControlBlock = repmat(sparse(diag(ones(N,1))),Nu,Nu);
    
    TimeStateBlock   = ones(size(StateStateBlock,1),2);
    TimeControlBlock = ones(size(StateControlBlock,1),2);
    StateTimeBlock   = TimeStateBlock.';
    ControlTimeBlock = TimeControlBlock.';
    
    rowcol = shift+1:shift+MAIN.nPhaseDecVar(iphase);
    
    Htotal(rowcol,rowcol) = [StateStateBlock,ControlStateBlock,TimeStateBlock;...
        StateControlBlock,ControlControlBlock,TimeControlBlock;
        StateTimeBlock,ControlTimeBlock,ones(2)];
    
    shift = rowcol(end,end);
end
Htotal = tril(Htotal);
end