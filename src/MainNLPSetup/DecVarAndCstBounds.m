function problem = DecVarAndCstBounds(problem)
%******************************************************************%
% This function determines the lower and upper bounds on all of the
% decision variables and all of the constraints.
%
% The decision variables are the time, state and control.
% The constraints are the dynamic defects, quadrature defects,
% path constraints, boundary constraints and linkages 
%
% Inputs:
% - problem :       User supplied problem struct plus 
%                   new fields from previous functions
%
% Outputs:
% - problem.DecVarLB    :   Lower Bounds on Decision Variables
% - problem.DecVarUB    :   Upper Bounds on Decision Variables
% - problem.CstLB       :   Lower Bounds on Constraints
% - problem.CstUB       :   Upper Bounds on Constraints
% - problem.nPhaseDecVar:   Vector of no. of Dec Vars per phase
% - problem.nPhaseCsts  :   Vector of no. of Constraints per phase
% - problem.nNonLinCsts :   Total number of Constraints
% - problem.nPhaseCstsLessBnd :  Vector of no. of Constraints minus 
%                               boundary constraints per phase
% - problem.nLinkCsts   :   Total number of Linkage Constraints
% - problem.DecVarIdx   :   Decision Variable Indexes in 
%                           Each Phase (cell)
% - problem.CstIdx      :   Constraint Indexes in Each Phase (cell)
% - problem.LinkCstIdx  :   Linkage Constraint Indexes
% - problem.Idx         :   State, Control, Time indexes as
%                           struct fields
% - problem.LinkJac.Jlinear : The jacobian of the linkage constraints.
%                             Since these are linear constraints, 
%                             the jacobian entries are +/- 1. 
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%******************************************************************%

nodes       = problem.nodes;
bounds      = problem.bounds.phase;
sizes       = problem.sizes;
nphases     = problem.nphases;
nStates     = sizes(1,:);
nControls   = sizes(2,:);
nPathCsts   = sizes(3,:);
nBndCsts    = sizes(4,:);
ShiftDecVar = 0;
DecVarCell  = cell(nphases,2);
CstCell     = cell(nphases,2);
nCell       = cell(nphases,3);
for iphase=1:nphases
    % Upper and Lower Bounds on initial and final time
    t0Low = bounds(iphase).initialtime.lb;
    t0Upp = bounds(iphase).initialtime.ub;
    tfLow = bounds(iphase).finaltime.lb;
    tfUpp = bounds(iphase).finaltime.ub;
    % Upper and lower bounds on all state cases (initial, during, final)
    if isequal(nStates(iphase),0)
        StateLBMatrix = [];
        StateUBMatrix = [];
    else
        [StateLBMatrix, StateUBMatrix] = ...
            BoundCollector(bounds,nodes,iphase,'state');
    end
    StateLBVector = StateLBMatrix(:);
    StateUBVector = StateUBMatrix(:);
    % Upper and lower bounds on control
    if isequal(nControls(iphase),0)
        ControlLBMatrix = [];
        ControlUBMatrix = [];
    else
        [ControlLBMatrix, ControlUBMatrix] = ...
            BoundCollector(bounds,nodes,iphase,'control');
    end
    ControlLBVector = ControlLBMatrix(:);
    ControlUBVector = ControlUBMatrix(:);
    % Upper and lower bounds on path constraints
    if isequal(nPathCsts,0)
        PathLBMatrix = [];
        PathUBMatrix = [];
    else
        [PathLBMatrix, PathUBMatrix] = ...
            BoundCollector(bounds,nodes,iphase,'path');
    end
    PathLBVector = PathLBMatrix(:);
    PathUBVector = PathUBMatrix(:);
    % Upper and lower bounds on boundary constraints
    if isequal(nBndCsts,0)
        BndLBVector = [];
        BndUBVector = [];
    else
        [BndLBVector, BndUBVector] = ...
            BoundCollector(bounds,nodes,iphase,'boundary');
    end
    % Pack up State, Control, Path Cst and Boundary Cst bounds into 
    % matrices and into column vectors. Also get the indicies for 
    % the decision variables (t,x,u) and the constraints
    
    % There are nodes*nStates amount of defect constraints
    DefectEqCstVector = zeros((nodes(iphase))*nStates(iphase),1);
    Z.phase(iphase).DecVar.lb   = [StateLBVector; ControlLBVector; t0Low; tfLow];
    Z.phase(iphase).DecVar.ub   = [StateUBVector; ControlUBVector; t0Upp; tfUpp];
    Z.phase(iphase).Cst.lb      = [DefectEqCstVector; PathLBVector; BndLBVector];
    Z.phase(iphase).Cst.ub      = [DefectEqCstVector; PathUBVector; BndUBVector];
    Z.phase(iphase).length.DecVar   = length(Z.phase(iphase).DecVar.lb);
    Z.phase(iphase).length.Cst      = length(Z.phase(iphase).Cst.lb);
    Z.phase(iphase).length.DefectPath =...
        length(Z.phase(iphase).Cst.lb) - length(BndLBVector);
    packDecVar.phase(iphase).stateIdx = (ShiftDecVar+1):...
        (ShiftDecVar+(nodes(iphase)+1)*nStates(iphase));  
    packDecVar.phase(iphase).controlIdx =...
        (packDecVar.phase(iphase).stateIdx(end)+1):...
        (packDecVar.phase(iphase).stateIdx(end)+...
        nodes(iphase)*nControls(iphase));
    
    if ~isempty(packDecVar.phase(iphase).controlIdx)
        t0Idx = packDecVar.phase(iphase).controlIdx(end)+1;
    else
        t0Idx = packDecVar.phase(iphase).stateIdx(end)+1;
    end
    
    tfIdx = t0Idx+1;
    packDecVar.phase(iphase).timeIdx = [t0Idx, tfIdx];
    ShiftDecVar             = ShiftDecVar+Z.phase(iphase).length.DecVar;
    DecVarCell{iphase,1}    = Z.phase(iphase).DecVar.lb;
    DecVarCell{iphase,2}    = Z.phase(iphase).DecVar.ub;
    CstCell{iphase,1}       = Z.phase(iphase).Cst.lb;
    CstCell{iphase,2}       = Z.phase(iphase).Cst.ub;
    nCell{iphase,1}         = Z.phase(iphase).length.DecVar;
    nCell{iphase,2}         = Z.phase(iphase).length.Cst;
    nCell{iphase,3}         = Z.phase(iphase).length.DefectPath;
end

% Get link constraint upper and lower bounds as column vectors
if nphases - 1 > 0 
    linkLBVector = [problem.bounds.link.lb].';
    linkUBVector = [problem.bounds.link.ub].';
else
    linkLBVector = [];
    linkUBVector = [];
end

nLinkCst = length(linkLBVector);
DecVarLB = vertcat(DecVarCell{:,1});
DecVarUB = vertcat(DecVarCell{:,2});
CstLBVector = vertcat(CstCell{:,1});
CstUBVector = vertcat(CstCell{:,2});
CstLB = [CstLBVector; linkLBVector];
CstUB = [CstUBVector; linkUBVector];

% Lower and Upper Bounds on Decision Variables in every phase
% These will be passed into ipopt options struct
problem.DecVarLB = DecVarLB;
problem.DecVarUB = DecVarUB;
% Lower and Upper Bounds on all Constraints in every phase
% These will be passed into ipopt options struct
problem.CstLB = CstLB;
problem.CstUB = CstUB;
% Number of Decision Variables In Each Phase (row vector)
problem.nPhaseDecVar = horzcat(nCell{:,1});
% Number of Csts (excl. linkages) In Each Phase (row vector)
problem.nPhaseCsts = horzcat(nCell{:,2});
% Number of Csts (excl. Bnd Csts and Linkages) In Each Phase (row vector)
problem.nPhaseCstsLessBnd = horzcat(nCell{:,3});
% Total number of Constraints (incl. linkages)
problem.nNonLinCsts = length(CstLB);
% Total number of Linkage Constraints
problem.nLinkCsts = nLinkCst;    
% Struct of Individual Decision Variable Indexes (State, Control, Time) 
problem.Idx = packDecVar;

% Compute the jacobian values for linkage constraints.
% Since these are all linear constraints, entries are +/- 1.
% Saves time in the actual jacobian calculation.
nlinks = nphases-1;
if nlinks > 0
    Jlinear = sparse(nLinkCst,sum(problem.nPhaseDecVar));
    Nx = nStates(1);
    NCsts = sum(problem.nPhaseCsts);
    problem.LinkJac.rowIdx = NCsts+1:NCsts+nlinks*(Nx+1);
    TimeRow = [1:Nx+1:(nlinks)*(Nx+1)];    
    Variables = [0,cumsum(problem.nPhaseDecVar)];
    for ilink = 1:nlinks
        N0 = nodes(ilink+1);
        Nf = nodes(ilink);
        StateRow = TimeRow(ilink)+[1:Nx];
        t0col  = Variables(ilink+2)-1;
        tfcol  = Variables(ilink+1);
        x0col  = Variables(ilink+1)+[1:(N0+1):(N0+1)*Nx];
        xfcol  = [(Nf+1):(Nf+1):(Nf+1)*Nx]+Variables(ilink);
        Jlinear(TimeRow(ilink),[tfcol,t0col])=[-1,1];
        Jlinear(StateRow,[xfcol,x0col]) = [-1*speye(Nx),speye(Nx)];
    end
    problem.LinkJac.Jlinear = Jlinear;
end
end


function [LBMatrix, UBMatrix] = BoundCollector(bounds,nodes,iphase,Case)

switch Case
    
    case{'state'}
        x0LB   = bounds(iphase).initialstate.lb;
        x0UB   = bounds(iphase).initialstate.ub;
        xLB    = bounds(iphase).state.lb;
        xUB    = bounds(iphase).state.ub;
        xfLB   = bounds(iphase).finalstate.lb;
        xfUB   = bounds(iphase).finalstate.ub;
        LBMatrix = [x0LB;xLB.*ones(nodes(iphase)-1,1);xfLB];
        UBMatrix = [x0UB;xUB.*ones(nodes(iphase)-1,1);xfUB];
        
    case{'control'}
        uLB = bounds(iphase).control.lb;
        uUB = bounds(iphase).control.ub;
        
        if ~isempty(uLB)
            LBMatrix = uLB.*ones(nodes(iphase),1);
        else
            LBMatrix = [];
        end
        
        if ~isempty(uUB)
            UBMatrix = uUB.*ones(nodes(iphase),1);
        else
            UBMatrix = [];
        end

    case{'path'}
        pathLB = bounds(iphase).path.lb;
        pathUB = bounds(iphase).path.ub;
        if ~isempty(pathLB)
            LBMatrix = pathLB.*ones(nodes(iphase),1);
        else
            LBMatrix = [];
        end
        
        if ~isempty(pathUB)
            UBMatrix = pathUB.*ones(nodes(iphase),1);
        else
            UBMatrix = [];
        end
        
    case{'boundary'}
        LBMatrix = bounds(iphase).boundary.lb.';
        UBMatrix = bounds(iphase).boundary.ub.';
        
end
end