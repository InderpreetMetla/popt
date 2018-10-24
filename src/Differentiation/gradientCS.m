function gtotal = gradientCS(x)
%*********************************************************************%
% This function uses Complex Step Differentiation to compute
% the gradient of the objective function
%
% Inputs: 
%    - x   : Vector of decision variables
%
% Outputs: 
%    - g   : gradient of boundary and path objectives
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%*********************************************************************%
global MAIN

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
% Struct of Indexes of the State (x), Control (u) and Time (t) 
Idx_xut = MAIN.Idx;
% GPM Collocation Information
LGColloc = MAIN.LGColloc; 
% Auxdata
auxdata = MAIN.auxdata; 
% Unscale
x = (x-MAIN.scaling.decvar_shifts)./MAIN.scaling.decvar_scales;
% Compute Pertubations
dh = MAIN.derivatives.first.stepsize.gradient;
h = dh*(1+abs(x));
% Assign short variables to optimal control functions
phi = MAIN.funcs.BndObj;
g   = MAIN.funcs.PathObj;
dMayer = [];
dLagrange = [];
for iphase = 1:nphases
    N = nodes(iphase);
    Nx = nStates(iphase);
    Nu = nControls(iphase);
    auxdata.iphase = iphase;
    t0 = x(Idx_xut.phase(iphase).timeIdx(1));
    tf = x(Idx_xut.phase(iphase).timeIdx(2));
    W = LGColloc.phase(iphase).Weights; % Weight already transposed
    ht0 = h(Idx_xut.phase(iphase).timeIdx(1));
    htf = h(Idx_xut.phase(iphase).timeIdx(2));
    Tau = [LGColloc.phase(iphase).Points;1];
    Time = 0.5 * (tf - t0) .* (Tau + 1) + t0;
    tgr	= Time(1:end-1); 
    htgr = dh*(1+abs(tgr));
    StateVector = x(Idx_xut.phase(iphase).stateIdx);
    StateMatrix = reshape(StateVector,N+1,Nx);
    hx = reshape(h(Idx_xut.phase(iphase).stateIdx),N+1,Nx);
    
    xgr = StateMatrix(1:end-1,:);
    x0 = StateMatrix(1,:);
    xf = StateMatrix(end,:);
    ControlVector = x(Idx_xut.phase(iphase).controlIdx);
    
    if nControls(iphase)>0
        ugr = reshape(ControlVector,N,Nu);
        hu = reshape(h(Idx_xut.phase(iphase).controlIdx),N,Nu);
    else
        ugr = [];
        hu = [];
    end
    
    %******************************************************************%
    % Derivative Function Calculations
    %******************************************************************%
    PObj = g(tgr,xgr,ugr,auxdata);
%     dLagrangeStates = [];
%     dMayerStates = [];
    dMayerStates = cell(1,Nx);
    dLagrangeStates = cell(1,Nx);
    ex = eye(Nx);
    for iState = 1:Nx
        % Mayer Cost wrt Initial and Final State
        dMayerStates{iState} = [imag(phi(t0,tf,x0+hx(1,:)*1i.*ex(iState,:),...
            xf,auxdata))./hx(1,iState),zeros(1,N-1),imag(phi(t0,tf,...
            x0,xf+hx(end,:)*1i.*ex(iState,:),auxdata))./hx(end,iState)];
        % Lagrange Cost wrt State
        dLagrangeStates{iState} = [0.5*(tf-t0)*(W.'.*imag(g(tgr,xgr+...
            hx(1:end-1,:)*1i.*ex(iState,:),ugr,auxdata))./...
            hx(1:end-1,iState)).',0];
    end
    
    dLagrangeControl = cell(1,Nu);
    eu = speye(Nu);
    for iControl = 1:Nu
        % Dynamics and Path Constraints wrt Control    
        dLagrangeControl{iControl} = 0.5*(tf-t0)*(W.'.*imag(g(tgr,xgr,...
            ugr + hu*1i.*eu(iControl,:),auxdata))./hu(:,iControl)).';
    end    
    
    dLdt  = imag(g(tgr + htgr.*1i,xgr,ugr,auxdata))./htgr;
    dMdt0 = imag(phi(t0 + ht0.*1i,tf,x0,xf,auxdata))./ht0;
    dMdtf = imag(phi(t0,tf+ htf.*1i,x0,xf,auxdata))./htf;
    
    % Form Gradient
    alpha = 0.5*(1-Tau(1:end-1));
    beta = 0.5*(1+Tau(1:end-1));
    dLdt0 = -0.5*W*PObj+0.5*(tf-t0)*W*(alpha.*dLdt);
    dLdtf = 0.5*W*PObj+0.5*(tf-t0)*W*(beta.*dLdt);
    dMayer = [dMayer,horzcat(dMayerStates{:}),zeros(1,N*Nu),dMdt0,dMdtf];
    dLagrange = [dLagrange,horzcat(dLagrangeStates{:}),...
        horzcat(dLagrangeControl{:}),dLdt0,dLdtf];
end
gtotal = sparse(dMayer+dLagrange);
gtotal = gtotal*MAIN.scaling.invV;
end