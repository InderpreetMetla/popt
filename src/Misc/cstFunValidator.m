function cstFunValidator(Z,auxdata)
%******************************************************************%
% This function checks that the user's Dynamics, Path Constraints,
% and Boundary Constratints functions are correctly set up.
%
% Inputs:
% - Z: A vector created from the user supplied guess.
% - auxdata: The auxdata from the user.
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%******************************************************************%

global MAIN;

if ~isfield(MAIN.funcs,'Dynamics') || isempty(MAIN.funcs.Dynamics)
    error('Error: Must supply a Dynamics function');
end

if ~(isa(MAIN.funcs.Dynamics,'char')) && ...
        ~(isa(MAIN.funcs.Dynamics,'function_handle'))
    error(strcat('Error: Invalid Dynamics function input.',...
        ' Must be a ''char'' or ''function_handle'''));
end

if ~isfield(MAIN.funcs,'PathCst') || isempty(MAIN.funcs.PathCst)
    error('Error: Must supply a path constraint function');
end

if ~(isa(MAIN.funcs.PathCst,'char')) && ...
        ~(isa(MAIN.funcs.PathCst,'function_handle'))
    error(strcat('Error: Invalid path constraint function input.',...
        ' Must be a ''char'' or ''function_handle'''));
end

if ~isfield(MAIN.funcs,'BndCst') || isempty(MAIN.funcs.BndCst)
    error('Error: Must supply a boundary constraint function');
end

if ~(isa(MAIN.funcs.BndCst,'char')) && ...
        ~(isa(MAIN.funcs.BndCst,'function_handle'))
    error(strcat('Error: Invalid boundary constraint function input.',...
        ' Must be a ''char'' or ''function_handle'''));
end

% Number of Nodes
Nodes = MAIN.nodes;
% Number of Phases
nphases = MAIN.nphases;
% Matrix of number of States, Controls, Path and Bnd Csts
sizes = MAIN.sizes;
% Number of States in Each Phase (row vector)
nStates = sizes(1,:);
% Number of Controls in Each Phase (row vector)
nControls = sizes(2,:);
% Struct of Indexes of the State (x), Control (u) and Time (t)
Idx_xut = MAIN.Idx;
% GPM Collocation Information
LGColloc = MAIN.LGColloc;
% Unscale
Z = (Z-MAIN.scaling.decvar_shifts)./MAIN.scaling.decvar_scales;
%******************************************************************%
for iphase = 1:nphases
    
    t0                  = Z(Idx_xut.phase(iphase).timeIdx(1));
    tf                  = Z(Idx_xut.phase(iphase).timeIdx(2));
    LGRPoints           = LGColloc.phase(iphase).Points;
    Tau                 = [LGRPoints;1];
    Time                = 0.5 * (tf - t0) .* (Tau + 1) + t0;
    TimeGR              = Time(1:end-1);
    
    StateVector         = Z(Idx_xut.phase(iphase).stateIdx);
    StateMatrix         = reshape(StateVector,...
        Nodes(iphase)+1,nStates(iphase));
    StateGRMatrix       = StateMatrix(1:end-1,:);
    x0                  = StateMatrix(1,:);
    xf                  = StateMatrix(end,:);
    
    ControlVector       = Z(Idx_xut.phase(iphase).controlIdx);
    
    if nControls(iphase)>0
        ControlMatrix = reshape(ControlVector,...
            Nodes(iphase),nControls(iphase));
    else
        ControlMatrix = [];
    end
    
    %******************************************************************%
    % Constraint Calculations
    %******************************************************************%
    auxdata.iphase = iphase;
    % Compute Dynamics
    DAE = MAIN.funcs.Dynamics(TimeGR,StateGRMatrix,ControlMatrix,auxdata);
    % Compute Path Constraints
    PathCst = MAIN.funcs.PathCst(TimeGR,StateGRMatrix,ControlMatrix,auxdata);
    % Compute Boundary Constraints
    BndCst = MAIN.funcs.BndCst(t0,tf,x0,xf,auxdata);
    
    %******************************************************************%
    % Construct Constraint Vector
    %******************************************************************%
    
    % Dynamics Input Validation
    [nRow,nCol] = size(DAE);
    if ~isempty(DAE)
        if ~isequal(nRow,Nodes(iphase))
            error(strcat('Error: Number of rows of the Dynamics',...
                ' Output in Phase %i is not correct. Should be',...
                ' #rows = #nodes.'),iphase);
        end
        if ~isequal(nCol,sizes(1,iphase))
            error(strcat('Error: Number of Dynamics Outputted',...
                ' in Phase %i Does Not Match',...
                ' the Number of Bounds Provided on States.',...
                ' Either check Provided Bounds or the Dynamics',...
                ' Function.'),iphase);
        end
    end
    
    % Path Constraint Calculation and Input Validation
    if ~isempty(PathCst)
        [nRow,nCol] = size(PathCst);
        if ~isequal(nRow,Nodes(iphase))
            error(strcat('Error: Number of rows of the Path',...
                ' Constraint Vector in Phase %i is not',...
                ' correct. Should be #rows = #nodes'),iphase);
        end
        if ~isequal(nCol,sizes(3,iphase))
            error(strcat('Error: Number of Path',...
                ' Constraints Outputted in Phase %i Does ',...
                ' Not Match the Number of Bounds Provided.',...
                ' Either check Provided Bounds or the Path',...
                ' Constraint Function.'),iphase);
        end
    end
    
    % Boundary Constraint Calculation and Input Validation
    if ~isempty(BndCst)
        [nRow, nCol] = size(BndCst);
        if ~isequal(nRow,1)
            error(strcat('Error: Boundary Constraint should be',...
                ' a single row vector ouput'));
        end
        if ~isequal(nCol,sizes(4,iphase))
            error(strcat('Error: Number of Boundary',...
                ' Constraints Outputted in Phase %i Does ',...
                ' Not Match the Number of Bounds Provided.',...
                ' Either check Provided Bounds or the Boundary',...
                ' Constraint Function.'),iphase);
        end
    end
    
end

end