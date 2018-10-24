function H = hessianFDScaled(Z,sigma,lambda)
%*************************************************************%
% This function uses Finite Differences to compute
% the Hessian of a Scaled NLP. It does so in a sparse manner by 
% actually deriving the optimal control functions rather than
% objFun and cstFun. This is more efficient. 
%
% Inputs: 
%    - Z       : Vector of decision variables
%    - sigma   : Scalar Lagrangian Hessian multiplier for the  
%                functions dependent only on boundary time and state. 
%    - lambda  : Vector Lagrangian Hessian multiplier for the 
%                functions dependent time, state and control.
% Outputs: 
%    - H       : Lagrangian Hessian for the Scaled NLP
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%*************************************************************%
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
% Number of Path Csts in Each Phase (row vector)
nPathCsts = sizes(3,:);
% Number of Bnd Csts in Each Phase (row vector)
nBndCsts = sizes(4,:);
% Struct of Indexes of the State (x), Control (u) and Time (t)
Idx_xut = MAIN.Idx;
% GPM Collocation Information
LGColloc = MAIN.LGColloc;
% Auxdata Struct
auxdata = MAIN.auxdata;
% Cum Sum the Number of Constraints (Excl. Linkages) In Each Phase
CstSum = cumsum(MAIN.nPhaseCsts);
% Number of Dec Var  In Each Phase
DecVarSum = MAIN.nPhaseDecVar;
rowshiftPhase = 0;
% Unscale for time multipliers
Zt = (Z-MAIN.scaling.decvar_shifts)./MAIN.scaling.decvar_scales;
% Compute Pertubations
dh = MAIN.derivatives.second.stepsize;
h = dh*(1+abs(Z));
% Assign short variables to optimal control functions
g = MAIN.funcs.PathObj;
a = MAIN.funcs.Dynamics;
c = MAIN.funcs.PathCst;
H = sparse(sum(MAIN.nPhaseDecVar),sum(MAIN.nPhaseDecVar));
HessShift = 0;
for iphase = 1:nphases
    N = nodes(iphase);
    Nx = nStates(iphase);
    Nu = nControls(iphase);
    Np = nPathCsts(iphase);
    Nb = nBndCsts(iphase);
    
    t0      = Z(Idx_xut.phase(iphase).timeIdx(1));
    tf      = Z(Idx_xut.phase(iphase).timeIdx(2));
    trange  = Zt(Idx_xut.phase(iphase).timeIdx(2))-...
        Zt(Idx_xut.phase(iphase).timeIdx(1));
    ht0     = h(Idx_xut.phase(iphase).timeIdx(1));
    htf     = h(Idx_xut.phase(iphase).timeIdx(2));
    Tau     = LGColloc.phase(iphase).Points;
    Time	= 0.5 * (tf - t0) .* ([Tau;1] + 1) + t0;
    alpha   = 0.5*(1-Tau);
    beta    = 0.5*(1+Tau);
    StateVector	= Z(Idx_xut.phase(iphase).stateIdx);
    StateMatrix	= reshape(StateVector,N+1,Nx);
    hx = reshape(h(Idx_xut.phase(iphase).stateIdx),N+1,Nx);
    hgr = hx(1:end-1,:);
    
    ControlVector = Z(Idx_xut.phase(iphase).controlIdx);
    if nControls(iphase)>0
        ControlMatrix = reshape(ControlVector,N,Nu);
        hu = reshape(h(Idx_xut.phase(iphase).controlIdx),N,Nu);
    else
        ControlMatrix = [];
        hu = [];
    end
    
    %**************************************************************%
    % Setup Inputs
    %**************************************************************%
    auxdata.iphase = iphase;
    w = LGColloc.phase(iphase).Weights;
    
    %**************************************************************%
    % Lambda Separation into gamma, phi and nu
    %**************************************************************%
    row = rowshiftPhase+1:rowshiftPhase+N*Nx;
    gamma = lambda(row).';
    row = rowshiftPhase+N*Nx+1:rowshiftPhase+N*(Nx+Np);
    phi = lambda(row).';
    row = rowshiftPhase+N*(Nx+Np)+1:rowshiftPhase+N*(Nx+Np)+Nb;
    nu = lambda(row).';
    rowshiftPhase = CstSum(iphase);
    
    %**************************************************************%
    % Pre-Allocations
    %**************************************************************%
    dLExx = cell(Nx);
    dLEt0x = cell(1,Nx);
    dLEtfx = cell(1,Nx);
    dLIxx = cell(Nx);
    dLIxu = cell(Nu,Nx);
    dLIuu = cell(Nu);
    dLIt0x = cell(1,Nx);
    dLItfx = cell(1,Nx);
    dLIt0u = cell(1,Nu);
    dLItfu = cell(1,Nu);

    %**************************************************************%
    % Create Separate Lagrangians
    %**************************************************************%
    t = Time;
    tgr = t(1:end-1);
    htgr = dh*(1+abs(tgr));
    x = StateMatrix;
    xgr = x(1:end-1,:);
    u = ControlMatrix;
    x0 = x(1,:); xf = x(end,:);
    LE0 = LE(t0,tf,x0,xf,auxdata,sigma,nu);
    % Path Objective = po
    po = PorC(tgr,xgr,u,auxdata,g,N);
    % Path Constraint = pc
    pc = PorC(tgr,xgr,u,auxdata,c,N);
    % Dynamics = dyn
    dyn = d(tgr,xgr,u,auxdata,a,N);
    
    %**************************************************************%
    % Compute Derivatives
    %**************************************************************%
    LEt02p = LE(t0+2.00*ht0,tf,x0,xf,auxdata,sigma,nu);
    LEtf2p = LE(t0,tf+2.00*htf,x0,xf,auxdata,sigma,nu);
    LEt0tf = LE(t0+ht0,tf+htf,x0,xf,auxdata,sigma,nu);
    LEt0m = LE(t0+ht0,tf,x0,xf,auxdata,sigma,nu);
    LEtfm = LE(t0,tf+htf,x0,xf,auxdata,sigma,nu);
    % Partials of LE wrt t0 and tf
    dLEt0t0 = (LEt02p-2.00*LEt0m+LE0)./ht0.^2;
    dLEtftf = (LEtf2p-2.00*LEtfm+LE0)./htf.^2;
    dLEt0tf = (LEt0tf-LEt0m-LEtfm+LE0)./(ht0.*htf);
    
    ex = speye(Nx);
    eu = speye(Nu);
    
    for i = 1:Nx
        for j = 1:i
            LEx0ijp = LE(t0,tf,x0+hx(1,:).*(ex(i,:)+ex(j,:)),xf,...
                auxdata,sigma,nu);
            LExfijp = LE(t0,tf,x0,xf+hx(end,:).*(ex(i,:)+ex(j,:)),...
                auxdata,sigma,nu);
            LEx0im = LE(t0,tf,x0+hx(1,:).*ex(i,:),xf,...
                auxdata,sigma,nu);
            LExfim = LE(t0,tf,x0,xf+hx(end,:).*ex(i,:),...
                auxdata,sigma,nu);
            LEx0jm = LE(t0,tf,x0+hx(1,:).*ex(j,:),xf,...
                auxdata,sigma,nu);
            LExfjm = LE(t0,tf,x0,xf+hx(end,:).*ex(j,:),...
                auxdata,sigma,nu);
            LExfix0j = LE(t0,tf,x0+hx(1,:).*ex(j,:),...
                xf+hx(end,:).*+ex(i,:),auxdata,sigma,nu);
            
            % Partials of LE wrt x0 and xf
            dLEx0ix0j = (LEx0ijp-LEx0im-LEx0jm+LE0)./...
                (hx(1,i).*hx(1,j));
            dLExfixfj = (LExfijp-LExfim-LExfjm+LE0)./...
                (hx(end,i).*hx(end,j));
            dLExfix0j = (LExfix0j-LEx0jm-LExfim+LE0)./...
                (hx(1,j).*hx(end,i));
            daxixj = (d(tgr,xgr+hgr.*ex(i,:)+hgr.*ex(j,:),u,auxdata,a,N)-...
                d(tgr,xgr+hgr.*ex(i,:),u,auxdata,a,N)-...
                d(tgr,xgr+hgr.*ex(j,:),u,auxdata,a,N)+dyn)./(hgr(:,i).*hgr(:,j));
            if ~isempty(po)
                dgxixj = (PorC(tgr,xgr+hgr.*ex(i,:)+hgr.*ex(j,:),u,auxdata,g,N)-...
                    PorC(tgr,xgr+hgr.*ex(i,:),u,auxdata,g,N)-...
                    PorC(tgr,xgr+hgr.*ex(j,:),u,auxdata,g,N)+po)./(hgr(:,i).*hgr(:,j));
            else
                dgxixj = zeros(N,1);
            end
            if ~isempty(pc)
                dcxixj = (PorC(tgr,xgr+hgr.*ex(i,:)+hgr.*ex(j,:),u,auxdata,c,N)-...
                    PorC(tgr,xgr+hgr.*ex(i,:),u,auxdata,c,N)-...
                    PorC(tgr,xgr+hgr.*ex(j,:),u,auxdata,c,N)+pc)./(hgr(:,i).*hgr(:,j));
            else
                dcxixj = zeros(N,1);
            end
            % Partials of LE wrt xx
            dLExx{i,j} = [dLEx0ix0j,sparse(1,N);sparse(N-1,N+1);
                dLExfix0j,sparse(1,N-1),dLExfixfj];
            % Partials of LI wrt xx
            dLIxx{i,j} = [diag(0.5*trange*(sigma*w.'.*dgxixj - ...
                sum(reshape(gamma,N,Nx).*daxixj,2)) + sum(reshape(phi,...
                N,Np).*dcxixj,2)),sparse(N,1);sparse(1,N),0];
            if i~=j
                dLExx{j,i} = sparse(N+1,N+1);
                dLIxx{j,i} = sparse(N+1,N+1);
            end
        end
        
        LEt0x0i = LE(t0+ht0,tf,x0+hx(1,:).*ex(i,:),xf,auxdata,sigma,nu);
        LEtfx0i = LE(t0,tf+htf,x0+hx(1,:).*ex(i,:),xf,auxdata,sigma,nu);
        LEt0xfi = LE(t0+ht0,tf,x0,xf+hx(end,:).*ex(i,:),auxdata,sigma,nu);
        LEtfxfi = LE(t0,tf+htf,x0,xf+hx(end,:).*ex(i,:),auxdata,sigma,nu);
        LEx0m = LE(t0,tf,x0+hx(1,:).*ex(i,:),xf,auxdata,sigma,nu);
        LExfm = LE(t0,tf,x0,xf+hx(end,:).*ex(i,:),auxdata,sigma,nu);
        % Partials of LE wrt t0 or tf and x0 or xf
        dLEt0x0i = (LEt0x0i-LEt0m-LEx0m+LE0)./(ht0.*hx(1,i));
        dLEt0xfi = (LEt0xfi-LEt0m-LExfm+LE0)./(ht0.*hx(end,i));
        dLEtfx0i = (LEtfx0i-LEtfm-LEx0m+LE0)./(htf.*hx(1,i));
        dLEtfxfi = (LEtfxfi-LEtfm-LExfm+LE0)./(htf.*hx(end,i));
        % Partials of LE wrt t0 or x
        dLEt0x{1,i} = [dLEt0x0i, sparse(1,N-1), dLEt0xfi];
        dLEtfx{1,i} = [dLEtfx0i, sparse(1,N-1), dLEtfxfi];
        
        % Derive wrt Control
        for jj = 1:Nu
            daxiuj = (d(tgr,xgr+hgr.*ex(i,:),u+hu.*eu(jj,:),auxdata,a,N)-...
                d(tgr,xgr+hgr.*ex(i,:),u,auxdata,a,N) - d(tgr,xgr,...
                u+hu.*eu(jj,:),auxdata,a,N)+dyn)./(hgr(:,i).*hu(:,jj));
            if ~isempty(po)
                dgxiuj = (PorC(tgr,xgr+hgr.*ex(i,:),u+hu.*eu(jj,:),auxdata,g,N)-...
                    PorC(tgr,xgr+hgr.*ex(i,:),u,auxdata,g,N)-PorC(tgr,xgr,u+...
                    hu.*eu(jj,:),auxdata,g,N)+po)./(hgr(:,i).*hu(:,jj));
            else
                dgxiuj = zeros(N,1);
            end
            if ~isempty(pc)
                dcxiuj = (PorC(tgr,xgr+hgr.*ex(i,:),u+hu.*eu(jj,:),auxdata,c,N)-...
                    PorC(tgr,xgr+hgr.*ex(i,:),u,auxdata,c,N)-PorC(tgr,xgr,u+...
                    hu.*eu(jj,:),auxdata,c,N)+pc)./(hgr(:,i).*hu(:,jj));
            else
                dcxiuj = zeros(N,1);
            end
            % Partials of LI wrt xu
            dLIxu{jj,i} = [diag(0.5*trange*(sigma*w.'.*dgxiuj - ...
                sum(reshape(gamma,N,Nx).*daxiuj,2)) + sum(reshape(phi,...
                N,Np).*dcxiuj,2)),sparse(N,1)];
        end
        
        % Derive wrt Time
        dadx = (d(tgr,xgr+hgr.*ex(i,:),u,auxdata,a,N)-...
            d(tgr,xgr-hgr.*ex(i,:),u,auxdata,a,N))./(2*hgr(:,i));
        dgdx = (PorC(tgr,xgr+hgr.*ex(i,:),u,auxdata,g,N)-...
            PorC(tgr,xgr-hgr.*ex(i,:),u,auxdata,g,N))./(2*hgr(:,i));
        datxi = (d(tgr+htgr,xgr+hgr.*ex(i,:),u,auxdata,a,N)-...
            d(tgr,xgr+hgr.*ex(i,:),u,auxdata,a,N) - d(tgr+htgr,xgr,...
            u,auxdata,a,N)+dyn)./(hgr(:,i).*htgr);
        if ~isempty(po)
            dgtxi = (PorC(tgr+htgr,xgr+hgr.*ex(i,:),u,auxdata,g,N)-...
                PorC(tgr,xgr+hgr.*ex(i,:),u,auxdata,g,N) - PorC(tgr+htgr,xgr,...
                u,auxdata,g,N)+po)./(hgr(:,i).*htgr);
        else
            dgtxi = zeros(N,1);
        end
        if ~isempty(pc)
            dctxi = (PorC(tgr+htgr,xgr+hgr.*ex(i,:),u,auxdata,c,N)-...
                PorC(tgr,xgr+hgr.*ex(i,:),u,auxdata,c,N) - PorC(tgr+htgr,xgr,...
                u,auxdata,c,N)+pc)./(hgr(:,i).*htgr);
        else
            dctxi = zeros(N,1);
        end
        dLIt0x{1,i} = [(0.5*(sum(reshape(gamma,N,Nx).*dadx,2)-sigma*w.'.*dgdx)+...
            0.5*trange*alpha.*(sigma*w.'.*dgtxi-sum(reshape(gamma,N,Nx).*datxi,2))...
            +alpha.*(sum(reshape(phi,N,Np).*dctxi,2))).',0].*MAIN.scaling.tScales{iphase};
        dLItfx{1,i} = [(0.5*(-sum(reshape(gamma,N,Nx).*dadx,2)+sigma*w.'.*dgdx)+...
            0.5*trange*beta.*(sigma*w.'.*dgtxi-sum(reshape(gamma,N,Nx).*datxi,2))...
            +beta.*(sum(reshape(phi,N,Np).*dctxi,2))).',0].*MAIN.scaling.tScales{iphase};
    end
    % Partials of LE wrt u and x
    dLEux = sparse(N*Nu,(N+1)*Nx);
    % Partials of LE wrt uu
    dLEuu = sparse(N*Nu,N*Nu);
    % Partials of LE wrt t0 and u
    dLEt0u = sparse(1,N*Nu);
    % Partials of LE wrt tf and u
    dLEtfu = sparse(1,N*Nu);
    
    % Partial of LI wrt uu
    for i = 1:Nu
        for j = 1:i
                dauiuj = (d(tgr,xgr,u+hu.*eu(i,:)+hu.*eu(j,:),auxdata,a,N)-...
                    d(tgr,xgr,u+hu.*eu(i,:),auxdata,a,N)-...
                    d(tgr,xgr,u+hu.*eu(j,:),auxdata,a,N)+dyn)./(hu(:,i).*hu(:,j));
                if ~isempty(po)
                    dguiuj = (PorC(tgr,xgr,u+hu.*eu(i,:)+hu.*eu(j,:),auxdata,g,N)-...
                        PorC(tgr,xgr,u+hu.*eu(i,:),auxdata,g,N)-...
                        PorC(tgr,xgr,u+hu.*eu(j,:),auxdata,g,N)+po)./(hu(:,i).*hu(:,j));
                else
                    dguiuj = zeros(N,1);
                end
                if ~isempty(pc)
                    dcuiuj = (PorC(tgr,xgr,u+hu.*eu(i,:)+hu.*eu(j,:),auxdata,c,N)-...
                        PorC(tgr,xgr,u+hu.*eu(i,:),auxdata,c,N)-...
                        PorC(tgr,xgr,u+hu.*eu(j,:),auxdata,c,N)+pc)./(hu(:,i).*hu(:,j));
                else
                    dcuiuj = zeros(N,1);
                end
            dLIuu{i,j} = diag(0.5*trange*(sigma*w.'.*dguiuj - ...
                sum(reshape(gamma,N,Nx).*dauiuj,2)) + sum(reshape(phi,...
                N,Np).*dcuiuj,2));
            if j ~= i
                dLIuu{j,i} = sparse(N,N);
            end
        end
        % Derive wrt Time
        dadu = (d(tgr,xgr,u+hu.*eu(i,:),auxdata,a,N)-...
            d(tgr,xgr,u-hu.*eu(i,:),auxdata,a,N))./(2*hu(:,i));
        dgdu = (PorC(tgr,xgr,u+hu.*eu(i,:),auxdata,g,N)-...
            PorC(tgr,xgr,u-hu.*eu(i,:),auxdata,g,N))./(2*hu(:,i));
        datui = (d(tgr+htgr,xgr,u+hu.*eu(i,:),auxdata,a,N)-...
            d(tgr,xgr,u+hu.*eu(i,:),auxdata,a,N) - d(tgr+htgr,xgr,...
            u,auxdata,a,N)+dyn)./(hu(:,i).*htgr);
        if ~isempty(po)
            dgtui = (PorC(tgr+htgr,xgr,u+hu.*eu(i,:),auxdata,g,N)-...
                PorC(tgr,xgr,u+hu.*eu(i,:),auxdata,g,N) - PorC(tgr+htgr,xgr,...
                u,auxdata,g,N)+po)./(hu(:,i).*htgr);
        else
            dgtui = zeros(N,1);
        end
        if ~isempty(pc)
            dctui = (PorC(tgr+htgr,xgr,u+hu.*eu(i,:),auxdata,c,N)-...
                PorC(tgr,xgr,u+hu.*eu(i,:),auxdata,c,N) - PorC(tgr+htgr,xgr,...
                u,auxdata,c,N)+pc)./(hu(:,i).*htgr);
        else
            dctui = zeros(N,1);
        end
        dLIt0u{1,i} = ((0.5*(sum(reshape(gamma,N,Nx).*dadu,2)-sigma*w.'.*dgdu)+...
            0.5*trange*alpha.*(sigma*w.'.*dgtui-sum(reshape(gamma,N,Nx).*datui,2))...
            +alpha.*(sum(reshape(phi,N,Np).*dctui,2))).')*MAIN.scaling.tScales{iphase};
        dLItfu{1,i} = ((0.5*(-sum(reshape(gamma,N,Nx).*dadu,2)+sigma*w.'.*dgdu)+...
            0.5*trange*beta.*(sigma*w.'.*dgtui-sum(reshape(gamma,N,Nx).*datui,2))...
            +beta.*(sum(reshape(phi,N,Np).*dctui,2))).')*MAIN.scaling.tScales{iphase};
    end
    
    % LI wrt t0
    dadt = (d(tgr+htgr,xgr,u,auxdata,a,N)-d(tgr-htgr,xgr,u,auxdata,a,N))./(2*htgr);
    dgdt = (PorC(tgr+htgr,xgr,u,auxdata,g,N)-PorC(tgr-htgr,xgr,u,auxdata,g,N))./(2*htgr);
    datt = (d(tgr+2*htgr,xgr,u,auxdata,a,N)-2*d(tgr+htgr,xgr,u,auxdata,a,N)+...
        dyn)./(htgr.^2);
    if ~isempty(po)
        dgtt = (PorC(tgr+2*htgr,xgr,u,auxdata,g,N)-2*PorC(tgr+htgr,xgr,u,auxdata,g,N)+...
            po)./(htgr.^2);
    else
        dgtt = zeros(N,1);
    end
    if ~isempty(pc)
        dctt = (PorC(tgr+2*htgr,xgr,u,auxdata,c,N)-2*PorC(tgr+htgr,xgr,u,auxdata,c,N)+...
            pc)./(htgr.^2);
    else
        dctt = zeros(N,1);
    end
    
    dLIt0t0 = alpha.'*(sum(reshape(gamma,N,Nx).*dadt,2)-sigma*w.'.*dgdt)+...
        0.5*trange*alpha.'*((sigma*w.'.*dgtt-sum(reshape(gamma,N,...
        Nx).*datt,2)).*alpha)+alpha.'*((sum(reshape(phi,N,Np).*dctt,2)).*alpha);
    dLItft0 = 0.5*alpha.'*(-sum(reshape(gamma,N,Nx).*dadt,2)+sigma*w.'.*dgdt)+...
        0.5*beta.'*(sum(reshape(gamma,N,Nx).*dadt,2)-sigma*w.'.*dgdt)+...
        0.5*trange*alpha.'*((sigma*w.'.*dgtt-sum(reshape(gamma,N,...
        Nx).*datt,2)).*beta)+alpha.'*((sum(reshape(phi,N,Np).*dctt,2)).*beta);
    dLItftf = beta.'*(-sum(reshape(gamma,N,Nx).*dadt,2)+sigma*w.'.*dgdt)+...
        0.5*trange*beta.'*((sigma*w.'.*dgtt-sum(reshape(gamma,N,...
        Nx).*datt,2)).*beta)+beta.'*((sum(reshape(phi,N,Np).*dctt,2)).*beta);
    % Partials of LE
    dLE = [cell2mat(dLExx),sparse(Nx*(N+1),N*Nu),sparse(Nx*(N+1),1),...
        sparse(Nx*(N+1),1);dLEux,dLEuu,sparse(Nu*N,1),sparse(Nu*N,1);
            cell2mat(dLEt0x),dLEt0u,dLEt0t0,sparse(0);
            cell2mat(dLEtfx),dLEtfu,dLEt0tf,dLEtftf];
    % Partials of LE
    dLI = [cell2mat(dLIxx),sparse(Nx*(N+1),N*Nu), sparse(Nx*(N+1),1),sparse(Nx*(N+1),1);
        cell2mat(dLIxu),cell2mat(dLIuu),sparse(Nu*N,1),sparse(Nu*N,1);
        cell2mat(dLIt0x),cell2mat(dLIt0u),dLIt0t0*MAIN.scaling.tScales{iphase},...
        sparse(0);cell2mat(dLItfx),cell2mat(dLItfu),dLItft0*MAIN.scaling.tScales{iphase},...
        dLItftf*MAIN.scaling.tScales{iphase}];
    rowcolPhase = HessShift+1:HessShift+DecVarSum(iphase);
    H(rowcolPhase,rowcolPhase) = dLE + dLI;
    HessShift = rowcolPhase(end);
end
end

function output = LE(t0,tf,x0,xf,auxdata,sigma,nu)
global MAIN
iphase = auxdata.iphase;
tScales = MAIN.scaling.tScales{iphase};
tShifts = MAIN.scaling.tShifts{iphase};
xScales = MAIN.scaling.stateScales{iphase};
xShifts = MAIN.scaling.stateShifts{iphase};

t0 = (t0-tShifts).*tScales;
tf = (tf-tShifts).*tScales;
x0 = (x0-xShifts).*xScales;
xf = (xf-xShifts).*xScales;

MayerCost = MAIN.funcs.BndObj(t0,tf,x0,xf,auxdata);
BndConstraints = MAIN.funcs.BndCst(t0,tf,x0,xf,auxdata);
output = sigma*MayerCost  + nu*BndConstraints(:);
end

function output = d(t,x,u,auxdata,original_function,N)
% x and u will be perturbed values and they will be scaled
% first we need to unscale them
global MAIN
iphase = auxdata.iphase;

tScales = repmat(MAIN.scaling.tScales{iphase},N,1);
tShifts = repmat(MAIN.scaling.tShifts{iphase},N,1);
xScales = repmat(MAIN.scaling.stateScales{iphase},N,1);
xShifts = repmat(MAIN.scaling.stateShifts{iphase},N,1);
uScales = repmat(MAIN.scaling.controlScales{iphase},N,1);
uShifts = repmat(MAIN.scaling.controlShifts{iphase},N,1);
t = (t - tShifts).*tScales;
x = (x - xShifts).*xScales;
u = (u - uShifts).*uScales;
output = original_function(t,x,u,auxdata)./xScales;
end

function output = PorC(t,x,u,auxdata,original_function,N)
% x and u will be perturbed values and they will be scaled
% first we need to unscale them
global MAIN
iphase = auxdata.iphase;
tScales = repmat(MAIN.scaling.tScales{iphase},N,1);
tShifts = repmat(MAIN.scaling.tShifts{iphase},N,1);
xScales = repmat(MAIN.scaling.stateScales{iphase},N,1);
xShifts = repmat(MAIN.scaling.stateShifts{iphase},N,1);
uScales = repmat(MAIN.scaling.controlScales{iphase},N,1);
uShifts = repmat(MAIN.scaling.controlShifts{iphase},N,1);
t = (t - tShifts).*tScales;
x = (x - xShifts).*xScales;
u = (u - uShifts).*uScales;
output = original_function(t,x,u,auxdata);
end