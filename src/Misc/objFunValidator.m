function objFunValidator(Z,auxdata)
%******************************************************************%
% This function checks that the Path objective function and the
% Boundary objective function are correctly set up.
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

if ~isfield(MAIN.funcs,'PathObj') || isempty(MAIN.funcs.PathObj)
    error('Error: Must supply an path objective function');
end

if ~(isa(MAIN.funcs.PathObj,'char')) && ...
        ~(isa(MAIN.funcs.PathObj,'function_handle'))
    
    error(strcat('Error: Invalid objective function input.',...
        ' Must be a ''char'' or ''function_handle'''))
end

if ~isfield(MAIN.funcs,'BndObj') || isempty(MAIN.funcs.BndObj)
    error('Error: Must supply an boundary objective function');
end

if ~(isa(MAIN.funcs.BndObj,'char')) && ...
        ~(isa(MAIN.funcs.BndObj,'function_handle'))
    
    error(strcat('Error: Invalid objective function input.',...
        ' Must be a ''char'' or ''function_handle'''))
end

% Number of Nodes
Nodes	  = MAIN.nodes;
% Number of Phases
nphases	  = MAIN.nphases;
% Matrix of number of States, Controls, Path and Bnd Csts
sizes	  = MAIN.sizes;
% Number of States in Each Phase (row vector)
nStates	  = sizes(1,:);
% Number of Controls in Each Phase (row vector)
nControls = sizes(2,:);
% Struct of Indexes of the State (x), Control (u) and Time (t)
Idx_xut	  = MAIN.Idx;
% GPM Collocation Information
LGColloc  = MAIN.LGColloc;
% Unscale
Z = (Z-MAIN.scaling.decvar_shifts)./MAIN.scaling.decvar_scales;
%**************************************************************%
% Input Construction
%**************************************************************%
for iphase = 1:nphases
    
    t0 = Z(Idx_xut.phase(iphase).timeIdx(1));
    tf = Z(Idx_xut.phase(iphase).timeIdx(2));
    LGRPoints = LGColloc.phase(iphase).Points;
    Tau	= [LGRPoints;1];
    Time = 0.5 * (tf - t0) .* (Tau + 1) + t0;
    TimeGR = Time(1:end-1);
    
    StateVector = Z(Idx_xut.phase(iphase).stateIdx);
    StateMatrix = reshape(StateVector,...
        Nodes(iphase)+1,nStates(iphase));
    StateGRMatrix = StateMatrix(1:end-1,:);
    x0 = StateMatrix(1,:);
    xf = StateMatrix(end,:);
    
    ControlVector = Z(Idx_xut.phase(iphase).controlIdx);
    
    if nControls(iphase)>0
        ControlMatrix = reshape(ControlVector,...
            Nodes(iphase),nControls(iphase));
    else
        ControlMatrix = [];
    end
    
    %**************************************************************%
    % Cost Function Validation
    %**************************************************************%
    auxdata.iphase = iphase;
    % Compute Objective Cost. Output is a struct of Mayer
    % and Lagrange Cost
    
    PathObj = MAIN.funcs.PathObj(TimeGR,StateGRMatrix,ControlMatrix,auxdata);
    BndObj = MAIN.funcs.BndObj(t0,tf,x0,xf,auxdata);
    
    
    if isempty(BndObj)
        error(strcat('Error: No Boundary Objective specified in phase %i.',...
            ' It cannot be left empty. If there is no boundary cost,',...
            ' then in your boundary objective function set the',...
            ' output to ',...
            ' ''zeros(size(t0))''.',...
            ' Note: you must use this format to set the cost',...
            ' = 0 otherwise errors will occur.'), iphase);
    elseif isempty(PathObj)
        error(strcat('Error: No Path Objective specified in',...
            ' phase %i. It cannot be left empty. If there is no',...
            ' path cost, then in your path objective function',...
            ' set the output to',...
            ' ''zeros(size(t))''.',...
            ' Note: you must use this format to set',...
            ' the cost = 0 otherwise errors will',...
            ' occur.'), iphase);
    end
    
    if ~isscalar(BndObj)
        error('Error: Boundary Objective is not a scalar value in phase %i.',...
            iphase)
    end
    
    nRowCol = size(PathObj);
    if ~isequal(nRowCol(1), Nodes(iphase))
        error(strcat('Error: Path Objective in phase %i is not ',...
            ' a column vector with same number of rows as ',...
            ' there are nodes.'),iphase);
    elseif ~isequal(nRowCol(2),1)
        error(strcat('Error: Path Objective in phase %i is not ',...
            ' a single column vector.'),iphase);
    end
end

end