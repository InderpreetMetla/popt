function ConstraintVector = cstFun(Z)
%******************************************************************%
% This function computes the Constraints and is called by IPOPT. 
% Output is one column vector of constraint evaluations.
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%******************************************************************%

global MAIN;

% Number of Nodes 
nodes = MAIN.nodes;
% Number of Phases
nphases = MAIN.nphases;
% Matrix of number of States, Controls, Path and Bnd Csts
sizes = MAIN.sizes;
% Number of States in Each Phase (row vector)
nStates = sizes(1,:);
% Number of Controls in Each Phase (row vector)
nControls = sizes(2,:);
% Number of Pairs of Linkages
nlinks = nphases - 1;
% Struct of Indexes of the State (x), Control (u) and Time (t) 
Idx_xut = MAIN.Idx;
% GPM Collocation Information
LGColloc = MAIN.LGColloc; 
% Auxdata
auxdata = MAIN.auxdata;
%******************************************************************%
% Initialise the input for the cost calculation
ConstraintCell  = cell(nphases,1);
% Unscale
Z = (Z-MAIN.scaling.decvar_shifts)./MAIN.scaling.decvar_scales;
for iphase = 1:nphases
        
    t0 = Z(Idx_xut.phase(iphase).timeIdx(1));
    tf = Z(Idx_xut.phase(iphase).timeIdx(2));
    LGPoints = LGColloc.phase(iphase).Points;
    Tau  = [LGPoints;1];
    Time = 0.5 * (tf - t0) .* (Tau + 1) + t0;
    TimeGR = Time(1:end-1);
    
    StateVector = Z(Idx_xut.phase(iphase).stateIdx);
    StateMatrix = reshape(StateVector,nodes(iphase)+1,nStates(iphase));
    StateGRMatrix = StateMatrix(1:end-1,:);
    x0 = StateMatrix(1,:);
    xf = StateMatrix(end,:);
    
    ControlVector = Z(Idx_xut.phase(iphase).controlIdx);
    if nControls(iphase)>0
        ControlMatrix = reshape(ControlVector,...
            nodes(iphase),nControls(iphase));
    else
        ControlMatrix = [];
    end
    %**************************************************************%
    % Store Some Boundary Inputs for Linkage Constraints
    %**************************************************************%
    inputBnd.phase(iphase).initialtime = t0;
    inputBnd.phase(iphase).finaltime = tf;
    inputBnd.phase(iphase).initialstate = x0;
    inputBnd.phase(iphase).finalstate = xf;
    
    %**************************************************************%
    % Construct Constraint Vector
    %**************************************************************%
    auxdata.iphase = iphase;
    % Calculate Parameters associated with LG Collocation
    LGDiffMatrix = LGColloc.phase(iphase).DM;
    % LGDiffMatrix = LGColloc.phase(iphase).DM_Diag;
    % Dynamics Calculation
    F = MAIN.funcs.Dynamics(TimeGR,StateGRMatrix,ControlMatrix,auxdata);
    % Defects:
    Defects = LGDiffMatrix * StateMatrix - 0.5 * (tf - t0) * F;
%     Defects = Defects./MAIN.scaling.stateScales{iphase};
    % Path Constraint Calculation and Input Validation  
    PathConstraints = MAIN.funcs.PathCst(TimeGR,StateGRMatrix,...
        ControlMatrix,auxdata);
    % Boundary Constraint Calculation and Input Validation
    BndConstraints = MAIN.funcs.BndCst(t0,tf,x0,xf,auxdata);
    % Total Constraint Vectors (Appended into a Cell)
    ConstraintCell{iphase} = [Defects(:);PathConstraints(:);...
        BndConstraints(:)];
end
% Concatenate all constraints into single column vector
ConstraintVector = vertcat(ConstraintCell{:,1});

%******************************************************************%
% Construct Linkage Constraint Vector
%******************************************************************%
Linkages = [];
for ilink = 1:nlinks
    
    t0 = inputBnd.phase(ilink+1).initialtime;
    tf = inputBnd.phase(ilink).finaltime;
    TimeLinks = t0 - tf;
    
    x0 = inputBnd.phase(ilink+1).initialstate;
    xf = inputBnd.phase(ilink).finalstate;
    StateLinks = x0 - xf;
    
    Linkages = [Linkages,TimeLinks,StateLinks];
end

LinkageVector = Linkages(:);
ConstraintVector = [ConstraintVector; LinkageVector];

ConstraintVector = MAIN.scaling.W*ConstraintVector;
end