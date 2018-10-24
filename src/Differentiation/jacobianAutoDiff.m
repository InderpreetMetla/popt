function Jtotal = jacobianAutoDiff(x)
%*********************************************************************%
% This function uses the auto differentiation package MatlabAutoDiff
% by martinReasearch to compute the Jacobian of the NLP.
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
global MAIN TempMAIN
TempMAIN = MAIN;
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
% Number of Defect + Path Constraints In Each Phase (row vector)
nDefectPathPhaseCsts = MAIN.nPhaseCstsLessBnd;
% Number of Constraints 
nNonLinCsts = MAIN.nNonLinCsts;
% Struct of Indexes of the State (x), Control (u) and Time (t) 
Idx_xut = MAIN.Idx;
% GPM Collocation Information
LGColloc = MAIN.LGColloc; 
% Preallocate Jacobian
Jtotal = sparse(nNonLinCsts,sum(nPhaseDecVar));
% Unscale
x = (x-MAIN.scaling.decvar_shifts)./MAIN.scaling.decvar_scales;

rowshift = 0;
colshiftPhase = 0;

for iphase = 1:nphases
    N = nodes(iphase);
    Nx = nStates(iphase);
    Nu = nControls(iphase);
    Np = nPathCsts(iphase);
    Nb = nBndCsts(iphase);
    MAIN.N = N;
    MAIN.Nx = Nx;
    MAIN.Nu = Nu;
    MAIN.Np = Np;
    MAIN.Nb = Nb;
    MAIN.iphase = iphase; 
    t0 = x(Idx_xut.phase(iphase).timeIdx(1));
    tf = x(Idx_xut.phase(iphase).timeIdx(2));
    Tau = LGColloc.phase(iphase).Points;
    Time = 0.5*(tf - t0).*([Tau;1]+1)+t0;
    StateVector = x(Idx_xut.phase(iphase).stateIdx);
    StateMatrix = reshape(StateVector,N+1,Nx);
    ControlVector = x(Idx_xut.phase(iphase).controlIdx);
    if Nu > 0
        ControlMatrix = reshape(ControlVector,N,Nu);
    else
        ControlMatrix = [];
    end
    % DAE and Path Constraint Function Input Initialization
    MAIN.time     = Time(1:end-1);
    MAIN.state    = StateMatrix(1:end-1,:);
    MAIN.control  = ControlMatrix;
    % Boundary Constraint Function Input Initialization
    MAIN.initialtime = t0;
    MAIN.finaltime = tf;
    MAIN.initialstate = StateMatrix(1,:);
    MAIN.finalstate = StateMatrix(end,:);
    %******************************************************************%
    % Optimal Control Function Evaluation
    %******************************************************************%
    MAIN.auxdata.iphase = iphase;
    F = MAIN.funcs.Dynamics(Time(1:end-1),StateMatrix(1:end-1,:),...
        ControlMatrix,MAIN.auxdata);

    %******************************************************************%
    % Derivative Function Calculations
    %******************************************************************%
    TempMAIN = MAIN;
    % Dynamics and Path Constraints wrt State
    if Nx > 0 && MAIN.dependency(iphase).dynamics.state
        dFdx = AutoDiffJacobianAutoDiff(@DAEwrtX,StateVector);
    else
        dFdx = sparse(Nx*N,Nx*(N+1));
    end 
    if Np > 0 && MAIN.dependency(iphase).pathcst.state
        dPdx = AutoDiffJacobianAutoDiff(@PathwrtX,StateVector);
    else
        dPdx = sparse(Np*N,Nx*(N+1));
    end
    
    % Reset
    TempMAIN = MAIN;
    % Dynamics and Path Constraints wrt Control
    if Nu > 0 && MAIN.dependency(iphase).dynamics.control
        dFdu = AutoDiffJacobianAutoDiff(@DAEwrtU,ControlVector);
    else
        dFdu = sparse(Nx*N,Nu*N);
    end 
    TempMAIN = MAIN;
    if Np > 0 && MAIN.dependency(iphase).pathcst.control
        dPdu = AutoDiffJacobianAutoDiff(@PathwrtU,ControlVector);
    else
        dPdu = sparse(Np*N,Nu*N);
    end
    
    % Reset
    TempMAIN = MAIN;
    % Dynamics and Path Constraints wrt Time
    if MAIN.dependency(iphase).dynamics.time
        dFdt = AutoDiffJacobianAutoDiff(@DAEwrtT,Time(1:end-1));
        dFdt = spdiags(dFdt, -(0:N:N*(Nx-1)));
    else
        dFdt = sparse(N,Nx);
    end
    if Np > 0 
        if MAIN.dependency(iphase).pathcst.time
            dPdt = AutoDiffJacobianAutoDiff(@PathwrtT,Time(1:end-1));
            dPdt = spdiags(dPdt, -(0:N:N*(Np-1)));
            dPdt0 = 0.5*(1-Tau).*dPdt;
            dPdtf = 0.5*(1+Tau).*dPdt;
        else
            dPdt0 = sparse(N*Np,1);
            dPdtf = sparse(N*Np,1);
        end
    else
        dPdt0 = [];
        dPdtf = [];
    end
  
    %************************************************************%
    % Insert Derivatives into Jacobian
    %************************************************************%
    DM = LGColloc.phase(iphase).DM;
    colshift = colshiftPhase;
    row = rowshift+1:rowshift+nDefectPathPhaseCsts(iphase);
    col = colshift+1:colshift+nPhaseDecVar(iphase);
    
    dFdt0 = 0.5*(F-0.5*(1-Tau)*(tf-t0).*dFdt);
    dFdtf = -0.5*(F+0.5*(1+Tau)*(tf-t0).*dFdt);
    
    Jtotal(row,col) = [kron(eye(Nx),DM)-0.5*(tf-t0)*dFdx,...
        -0.5*(tf-t0)*dFdu,dFdt0(:),dFdtf(:);...
        dPdx,dPdu,dPdt0(:),dPdtf(:)];
    
    rowshift = row(end);
    

    % Boundary Constraints
    if Nb > 0
        
        % Bnd Cst wrt t0
        if MAIN.dependency(iphase).bndcst.initialtime
            % Reset
            TempMAIN = MAIN;
            dBdt0 = AutoDiffJacobianAutoDiff(@Bndwrtt0,t0);
        else
            dBdt0 = sparse(Nb,1);
        end
        
        % Bnd Cst wrt tf
        if MAIN.dependency(iphase).bndcst.finaltime
            TempMAIN = MAIN;
            dBdtf = AutoDiffJacobianAutoDiff(@Bndwrttf,tf);
        else
            dBdtf = sparse(Nb,1);
        end
        % Bnd Cst wrt x0
        if MAIN.dependency(iphase).bndcst.initialstate
            TempMAIN = MAIN;
            dBdx0 = AutoDiffJacobianAutoDiff(@Bndwrtx0,StateMatrix(1,:));
        else
            dBdx0 = sparse(Nb,Nx);
        end
        % Bnd Cst wrt tf
        if MAIN.dependency(iphase).bndcst.finalstate
            TempMAIN = MAIN;
            dBdxf = AutoDiffJacobianAutoDiff(@Bndwrtxf,StateMatrix(end,:));
        else
            dBdxf = sparse(Nb,Nx);
        end
        row = rowshift+1:rowshift+Nb;
        col = colshiftPhase+...
            [(1:Nx)*(N+1)-N,(1:Nx)*(N+1),N*(Nx+Nu)+Nx+[1:2]];
        Jtotal(row,col) = [dBdx0,dBdxf,dBdt0,dBdtf];
    end
    rowshift = row(end);
    colshiftPhase = col(end);      

end
if nlinks > 0
    Jtotal(MAIN.LinkJac.rowIdx,:) = MAIN.LinkJac.Jlinear;
end
Jtotal = MAIN.scaling.W*Jtotal*MAIN.scaling.invV;
end

function Dynamics = DAEwrtX(Z)
global TempMAIN
% Number of Nodes 
N = TempMAIN.N;
% Number of States
Nx = TempMAIN.Nx;
StateMatrix = reshape(Z,N+1,Nx);
TempMAIN.state = StateMatrix(1:end-1,:);
Dynamics = TempMAIN.funcs.Dynamics(TempMAIN.time,TempMAIN.state,...
    TempMAIN.control,TempMAIN.auxdata);
end

function Dynamics = DAEwrtU(Z)
global TempMAIN
% Number of Nodes 
N = TempMAIN.N;
% Number of Controls
Nu = TempMAIN.Nu;
ControlMatrix = reshape(Z,N,Nu);
TempMAIN.control = ControlMatrix;
Dynamics = TempMAIN.funcs.Dynamics(TempMAIN.time,TempMAIN.state,...
    TempMAIN.control,TempMAIN.auxdata);
end

function Dynamics = DAEwrtT(Z)
global TempMAIN
% Time = 0.5*(Z(2)-Z(1)).*...
%     ([TempMAIN.LGColloc.phase(iphase).Points;1]+1)+Z(1); 
TempMAIN.time = Z;
Dynamics = TempMAIN.funcs.Dynamics(TempMAIN.time,TempMAIN.state,...
    TempMAIN.control,TempMAIN.auxdata);
end

function PathCst = PathwrtX(Z)
global TempMAIN
% Number of Nodes 
N = TempMAIN.N;
% Number of States
Nx = TempMAIN.Nx;
StateMatrix = reshape(Z,N+1,Nx);
TempMAIN.state = StateMatrix(1:end-1,:);
PathCst = TempMAIN.funcs.PathCst(TempMAIN.time,TempMAIN.state,...
    TempMAIN.control,TempMAIN.auxdata);
end

function PathCst = PathwrtU(Z)
global TempMAIN
% Number of Nodes 
N = TempMAIN.N;
% Number of Controls
Nu = TempMAIN.Nu;
ControlMatrix = reshape(Z,N,Nu);
TempMAIN.control = ControlMatrix;
PathCst = TempMAIN.funcs.PathCst(TempMAIN.time,TempMAIN.state,...
    TempMAIN.control,TempMAIN.auxdata);
end

function PathCst = PathwrtT(Z)
global TempMAIN
TempMAIN.time = Z;
PathCst = TempMAIN.funcs.PathCst(TempMAIN.time,TempMAIN.state,...
    TempMAIN.control,TempMAIN.auxdata);
end

function BndCst = Bndwrtt0(Z)
global TempMAIN
TempMAIN.initialtime = Z;
BndCst = TempMAIN.funcs.BndCst(TempMAIN.initialtime,TempMAIN.finaltime,...
    TempMAIN.initialstate,TempMAIN.finalstate,TempMAIN.auxdata);
end

function BndCst = Bndwrttf(Z)
global TempMAIN
TempMAIN.finaltime = Z;
BndCst = TempMAIN.funcs.BndCst(TempMAIN.initialtime,TempMAIN.finaltime,...
    TempMAIN.initialstate,TempMAIN.finalstate,TempMAIN.auxdata);
end

function BndCst = Bndwrtx0(Z)
global TempMAIN
TempMAIN.initialstate = Z;
BndCst = TempMAIN.funcs.BndCst(TempMAIN.initialtime,TempMAIN.finaltime,...
    TempMAIN.initialstate,TempMAIN.finalstate,TempMAIN.auxdata);
end

function BndCst = Bndwrtxf(Z)
global TempMAIN
TempMAIN.finalstate = Z;
BndCst = TempMAIN.funcs.BndCst(TempMAIN.initialtime,TempMAIN.finaltime,...
    TempMAIN.initialstate,TempMAIN.finalstate,TempMAIN.auxdata);
end