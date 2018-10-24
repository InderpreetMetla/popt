function problem = ScaleNlp(problem)
%**********************************************************************%
% This function scales the optimal control problem using the bounds
% provided by the user. All variables are scaled. The only function
% scaled is the dynamics function. Rest are "scaled" by 1.
%
% Inputs:
%
% - problem struct
%
% Outputs:
% - problem struct with:
%       - .scaling  :  Contains all scaling parameters for each
%                      variable and each function.
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%**********************************************************************%

nodes   = problem.nodes;
sizes   = problem.sizes;
bounds  = problem.bounds.phase;
nNonLin = problem.nNonLinCsts;
nphases = problem.nphases;
nStates = sizes(1,:);
nControls = sizes(2,:);
nPathCsts = sizes(3,:);
nBndCsts = sizes(4,:);

decvar_scales   = [];
decvar_shifts   = [];
function_scales = [];

stateScales   = cell(1,nphases);
stateShifts   = stateScales;
controlScales = cell(1,nphases);
controlShifts = controlScales;
tScales = cell(1,nphases);
tShifts = tScales;
tfScales = cell(1,nphases);
tfShifts = tfScales;

nlinks = nphases-1;

if ~isfield(problem,'autoscale')
    problem.autoscale = 'off';
end

if (isequal(problem.autoscale,'on') || isequal(problem.autoscale,1))
    for iphase=1:nphases
        %*******************************************%
        % Compute the Decision Variable Scales 
        % Using Maximum and Minimum Values of Time, 
        % State and Control Provided by User
        %********************************************%
        N = nodes(iphase);
        Nx = nStates(iphase);
        Nu = nControls(iphase);
        Np = nPathCsts(iphase);
        Nb = nBndCsts(iphase);
        
        %***************%
        % State
        %***************%
        xLow = bounds(iphase).state.lb;
        xUpp = bounds(iphase).state.ub;
        if ~isequal(xLow,xUpp)
            scale.x = (xUpp - xLow).*ones(N+1,Nx);
            shift.x = (0.5-xUpp./(xUpp-xLow)).*ones(N+1,1);
        else
            scale.x = ones(N+1,Nx);
            shift.x = zeros(N+1,1);
        end
        scale.x(scale.x==0) = 1;
        scale.x(abs(scale.x)==Inf) = ...
            sign(scale.x(abs(scale.x)==Inf))*1e10;
        shift.x(abs(shift.x)==Inf) = 0;
        shift.x(isnan(shift.x)) = 0;
        stateScales{iphase} = scale.x(1,:);
        stateShifts{iphase} = shift.x(1,:);
        
        %***************%
        % Control
        %***************%
        if Nu > 0
            uLow = bounds(iphase).control.lb;
            uUpp = bounds(iphase).control.ub;
        else
            uLow = [];
            uUpp = [];
        end
        if ~isequal(uLow,uUpp)
            scale.u = (uUpp-uLow).*ones(N,Nu);
            shift.u = (0.5-uUpp./(uUpp-uLow)).*ones(N,1);
        else
            scale.u = ones(N,Nu);
            shift.u = zeros(N,1);
        end
        scale.u(scale.u==0) = 1;
        scale.u(abs(scale.u)==Inf) = ...
            sign(scale.u(abs(scale.u)==Inf))*1e10;
        shift.u(abs(shift.u)==Inf) = 0;
        shift.u(isnan(shift.u)) = 0;
        controlScales{iphase} = scale.u(1,:);
        controlShifts{iphase} = shift.u(1,:);
        
        %***************%
        % Time
        %***************%
        t0Low = bounds(iphase).initialtime.lb;
        tfUpp = bounds(iphase).finaltime.ub;
        scale.t = tfUpp-t0Low;
        shift.t = 0.5-tfUpp/(tfUpp-t0Low);
        tScales{iphase} = scale.t;
        tShifts{iphase} = shift.t;
        
        %***************%
        % All Dec Vars
        %***************%
        decvar_scales = [decvar_scales;[scale.x(:);scale.u(:);...
            scale.t;scale.t]];
        decvar_shifts = [decvar_shifts;[shift.x(:);shift.u(:);...
            shift.t;shift.t]];

        %*******************************************%
        % Compute the Function Scales
        %*******************************************%
        dynmatrix = scale.x(1:end-1,:);
        scale.p = ones(N,Np);
        scale.b = ones(1,Nb);        
        function_scales = [function_scales;dynmatrix(:);...
            scale.p(:);scale.b(:)];
    end
    %*******************************************%
    % Compute the Linkage Scales
    %*******************************************%
    link_scales = [];
    for ilink = 1:nlinks
        Nx = nStates(1);
        scale.l = ones(1,Nx+1);
        link_scales = [link_scales;scale.l(:)];
    end
    function_scales = [function_scales;link_scales(:)];
else
    %*******************************************%
    % If no scaling required
    %*******************************************%
    for iphase=1:nphases
        Nx = nStates(iphase);
        Nu = nControls(iphase);
        stateScales{iphase} = ones(1,Nx);
        controlScales{iphase} = ones(1,Nu);
        stateShifts{iphase} = zeros(1,Nx);
        controlShifts{iphase} = zeros(1,Nu);
        tScales{iphase} = 1;
        tfScales{iphase} = 1;
        tShifts{iphase} = 0;
        tfShifts{iphase} = 0;
    end
    decvar_scales = ones(sum(problem.nPhaseDecVar),1);
    decvar_shifts = zeros(sum(problem.nPhaseDecVar),1);
    function_scales = ones(nNonLin,1);
end
%*******************************************%
% Append to problem
%*******************************************%
problem.scaling.stateScales     = stateScales;
problem.scaling.stateShifts     = stateShifts;

problem.scaling.controlScales   = controlScales;
problem.scaling.controlShifts   = controlShifts;

problem.scaling.tScales = tScales;
problem.scaling.tShifts = tShifts;

problem.scaling.decvar_scales   = 1./decvar_scales;
problem.scaling.decvar_shifts   = decvar_shifts;

% Compute a sparse diagonal matrix of the decision variable and 
% function scales - only need inverse of decision variable scales. 
% The inverse of a diagonal matrix is simply the reciprical of the
% diagonal elements.
problem.scaling.invV = spdiags(decvar_scales,0,...
    sum(problem.nPhaseDecVar),sum(problem.nPhaseDecVar));
problem.scaling.W = spdiags(1./function_scales,0,nNonLin,nNonLin);
end
