function Jtotal = jacobianCS(x)
%*********************************************************************%
% This function uses Complex Step Differentiation to compute the 
% Jacobian of the NLP.
% It does so in a sparse manner by deriving the optimal control
% functions rather than objFun and cstFun. This is more efficient. 
% 
% Inputs: 
%    - x   : Vector of decision variables
%
% Outputs: 
%    - Jtotal   : Jacobian of NLP
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
% Number of Path Constraints in Each Phase (row vector)
nPathCsts = sizes(3,:);
% Number of Boundary Constraints in Each Phase (row vector)
nBndCsts = sizes(4,:);
% Number of Pairs of Linkages
nlinks = nphases-1;
% Number of Decision Variables In Each Phase (row vector)
nPhaseDecVar = MAIN.nPhaseDecVar;
% Number of Constraints 
nNonLinCsts = MAIN.nNonLinCsts;
% Struct of Indexes of the State (x), Control (u) and Time (t) 
Idx_xut = MAIN.Idx;
% GPM Collocation Information
LGColloc = MAIN.LGColloc; 
% Auxdata
auxdata = MAIN.auxdata;
% Preallocate Jacobian
Jtotal = sparse(nNonLinCsts,sum(nPhaseDecVar));
% Unscale
x = (x-MAIN.scaling.decvar_shifts)./MAIN.scaling.decvar_scales;
% Compute Pertubations
dh = MAIN.derivatives.first.stepsize.jacobian;
h = dh*(1+abs(x));
% Row/Column Indicies
rowshift = 0;
colshiftPhase = 0;
% Assign short variables to optimal control functions
a = MAIN.funcs.Dynamics;
c = MAIN.funcs.PathCst;
b = MAIN.funcs.BndCst;
for iphase = 1:nphases
    N = nodes(iphase);
    Nx = nStates(iphase);
    Nu = nControls(iphase);
    Np = nPathCsts(iphase);
    Nb = nBndCsts(iphase);   
    t0 = x(Idx_xut.phase(iphase).timeIdx(1));
    tf = x(Idx_xut.phase(iphase).timeIdx(2));
    ht0 = h(Idx_xut.phase(iphase).timeIdx(1));
    htf = h(Idx_xut.phase(iphase).timeIdx(2));
    Tau = LGColloc.phase(iphase).Points;
    Time = 0.5 * (tf - t0) .* ([Tau;1] + 1) + t0;
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
        ControlMatrix = reshape(ControlVector,N,Nu);
        hu = reshape(h(Idx_xut.phase(iphase).controlIdx),N,Nu);
    else
        ControlMatrix = [];
        hu = [];
    end
    ugr = ControlMatrix;
    %******************************************************************%
    % Function Calculations
    %******************************************************************%
    auxdata.iphase = iphase;
    % Compute Dynamics
    F = a(tgr,xgr,ugr,auxdata);
    
    %******************************************************************%
    % Pre-Allocations
    %******************************************************************%
    dFdx = cell(1,Nx);
    dPdx = dFdx;
    dBdx0 = dFdx;
    dBdxf = dFdx;
    dFdu = cell(1,Nu);
    dPdu = dFdu;

    %******************************************************************%
    % Derivative Function Calculations
    %******************************************************************%
    ex = speye(Nx);
    % Dynamics and Path Constraints wrt State
    for iState = 1:Nx
        dFdx{iState} = imag(a(tgr,xgr+hx(1:end-1,:)*1i.*ex(iState,:),...
            ugr,auxdata))./hx(1:end-1,iState);
        if Np > 0 && MAIN.dependency(iphase).pathcst.state 
            dPdx{iState} = imag(c(tgr,xgr+hx(1:end-1,:)*1i.*ex(iState,:),...
                ugr,auxdata))./hx(1:end-1,iState);
        else
            dPdx{iState} = zeros(N,Np);
        end
        if Nb > 0
            if MAIN.dependency(iphase).bndcst.initialstate
                % Boundary Cst wrt x0
                dBdx0{iState} = imag(b(t0,tf,x0+hx(1,:)*1i.*ex(iState,:),...
                    xf,auxdata))./hx(1,iState);
            else
                dBdx0{iState} = sparse(1,Nb);
            end
            if MAIN.dependency(iphase).bndcst.finalstate
                dBdxf{iState} = imag(b(t0,tf,x0,xf+...
                    hx(end,:)*1i.*ex(iState,:),auxdata))./hx(end,iState);
            else
                dBdxf{iState} = sparse(1,Nb);
            end
        end
    end
    eu = speye(Nu);
    for iControl = 1:Nu
        % Dynamics and Path Constraints wrt Control
        dFdu{iControl} = imag(a(tgr,xgr,ugr + hu*1i.*eu(iControl,:),...
            auxdata))./hu(:,iControl);
        if Np > 0 && MAIN.dependency(iphase).pathcst.control
            dPdu{iControl} = ...
                imag(c(tgr,xgr,ugr + hu*1i.*eu(iControl,:),...
                auxdata))./hu(:,iControl);
        else
            dPdu{iControl} = zeros(N,Np);
        end
    end
    % Dynamics and Path Constraints wrt Time
    dFdt = imag(a(tgr+htgr.*1i,xgr,ugr,auxdata))./htgr;
    if Np > 0
        dPdt = imag(c(tgr+htgr.*1i,xgr,ugr,auxdata))./htgr;
    end

    if Nb > 0
        % Boundary Cst wrt t0
        if MAIN.dependency(iphase).bndcst.initialtime
            dBdt0 = imag(b(t0+ht0.*1i,tf,x0,xf,auxdata))./ht0;
        else
            dBdt0 = sparse(1,Nb);
        end
        % Boundary Cst wrt tf
        if MAIN.dependency(iphase).bndcst.finaltime
            dBdtf = imag(b(t0,tf+htf.*1i,x0,xf,auxdata))./htf;
        else
            dBdtf = sparse(1,Nb);
        end
    end
    
    %************************************************************%
    % Insert Derivatives into Jacobian
    %************************************************************%
    DM = LGColloc.phase(iphase).DM;
    col = colshiftPhase+1:colshiftPhase+Nx*(N+1)+N*Nu+2;
    for iState = 1:Nx
        Defect = [];
        row = rowshift+1:rowshift+N;
        for jState = 1:Nx
            Defect = [Defect,(iState == jState)*DM - ...
                0.5*(tf-t0)*[diag(dFdx{jState}(:,iState)),zeros(N,1)]];
        end
        for jControl = 1:Nu
            Defect = [Defect,-0.5*(tf-t0)*diag(dFdu{jControl}(:,iState))];
        end
        Defect =[Defect,0.5*(F(:,iState)-...
            (tf-t0)*0.5*(1-Tau).*dFdt(:,iState)),...
            -0.5*(F(:,iState)+(tf-t0)*0.5*(1+Tau).*dFdt(:,iState))];
        Jtotal(row,col) = sparse(Defect);
        rowshift = row(end);
    end
    
    for iPath = 1:Np
        row = rowshift+1:rowshift+N;
        Paths = [];
        for jState = 1:Nx
            Paths = [Paths,[diag(dPdx{jState}(:,iPath)),zeros(N,1)]];
        end
        for jControl = 1:Nu
            Paths = [Paths,diag(dPdu{jControl}(:,iPath))];
        end
        Paths = [Paths,0.5*(1-Tau).*dPdt(:,iPath),...
            0.5*(1+Tau).*dPdt(:,iPath)];
        Jtotal(row,col) = sparse(Paths);
        rowshift = row(end);
    end
    
    for iBnd = 1:Nb
        row = rowshift+iBnd;
%         Insert dB/dx0 and dB/dxf
        for jState = 1:Nx
            colx0 = colshiftPhase + jState*(N+1) - N;
            colxf = colshiftPhase + jState*(N+1);
            Jtotal(row,colx0) = sparse(dBdx0{jState}(iBnd));
            Jtotal(row,colxf) = sparse(dBdxf{jState}(iBnd));
        end
%         Insert dB/dt0 and dB/dtf
        col = colshiftPhase + Nx*(N+1)+Nu*N+[1:2];
        Jtotal(row,col) = sparse([dBdt0(iBnd),dBdtf(iBnd)]);
    end
    rowshift = row(end);
    colshiftPhase = col(end);
end    

% Linkage Constraints
if nlinks > 0
    Jtotal(MAIN.LinkJac.rowIdx,:) = MAIN.LinkJac.Jlinear;
end
Jtotal = MAIN.scaling.W*Jtotal*MAIN.scaling.invV;
end