function H = hessianFD(Z,sigma,lambda)
%*************************************************************%
% This function uses Finite Differences to compute
% the Hessian of the NLP. It does so in a sparse manner by 
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
%    - H       : Lagrangian Hessian
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
% Unscale
Z = (Z-MAIN.scaling.decvar_shifts)./MAIN.scaling.decvar_scales;
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
    ht0     = h(Idx_xut.phase(iphase).timeIdx(1));
    htf     = h(Idx_xut.phase(iphase).timeIdx(2));
    Tau     = LGColloc.phase(iphase).Points;
    Time	= 0.5 * (tf - t0) .* ([Tau;1] + 1) + t0;
    alpha   = 0.5*(1-Tau);
    beta    = 0.5*(1+Tau);
    %     hTime   = dh*(1+abs([Tau;1]));
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
    %     htgr = 2.05e-04;;
    x = StateMatrix;
    xgr = x(1:end-1,:);
    u = ControlMatrix;
    x0 = x(1,:); xf = x(end,:);
    LE0 = LE(t0,tf,x0,xf,auxdata,sigma,nu);
    % LI0 = LI(t,x,u,auxdata,sigma,gamma,phi,w,DM);
    % Path Objective = po
    po = g(tgr,xgr,u,auxdata);
    % Path Constraint = pc
    pc = c(tgr,xgr,u,auxdata);
    % Dynamics = dyn
    dyn = a(tgr,xgr,u,auxdata);
    
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
%             LEx0ixfj = LE(t0,tf,x0+hx(1,:).*ex(i,:),...
%                 xf+hx(end,:).*+ex(j,:),auxdata,sigma,nu);
            LExfix0j = LE(t0,tf,x0+hx(1,:).*ex(j,:),...
                xf+hx(end,:).*+ex(i,:),auxdata,sigma,nu);
            
            % Partials of LE wrt x0 and xf
            dLEx0ix0j = (LEx0ijp-LEx0im-LEx0jm+LE0)./...
                (hx(1,i).*hx(1,j));
            dLExfixfj = (LExfijp-LExfim-LExfjm+LE0)./...
                (hx(end,i).*hx(end,j));
%             dLEx0ixfj = (LEx0ixfj-LEx0im-LExfjm+LE0)./...
%                 (hx(1,i).*hx(end,j));
            dLExfix0j = (LExfix0j-LEx0jm-LExfim+LE0)./...
                (hx(1,j).*hx(end,i));
            
            daxixj = (a(tgr,xgr+hgr.*ex(i,:)+hgr.*ex(j,:),u,auxdata)-...
                a(tgr,xgr+hgr.*ex(i,:),u,auxdata)-...
                a(tgr,xgr+hgr.*ex(j,:),u,auxdata)+dyn)./(hgr(:,i).*hgr(:,j));
            if ~isempty(po)
                dgxixj = (g(tgr,xgr+hgr.*ex(i,:)+hgr.*ex(j,:),u,auxdata)-...
                    g(tgr,xgr+hgr.*ex(i,:),u,auxdata)-...
                    g(tgr,xgr+hgr.*ex(j,:),u,auxdata)+po)./(hgr(:,i).*hgr(:,j));
            else
                dgxixj = zeros(N,1);
            end
            if ~isempty(pc)
                dcxixj = (c(tgr,xgr+hgr.*ex(i,:)+hgr.*ex(j,:),u,auxdata)-...
                    c(tgr,xgr+hgr.*ex(i,:),u,auxdata)-...
                    c(tgr,xgr+hgr.*ex(j,:),u,auxdata)+pc)./(hgr(:,i).*hgr(:,j));
            else
                dcxixj = zeros(N,1);
            end
            % Partials of LE wrt xx
%             dLExx{i,j} = [dLEx0ix0j,sparse(1,N-1),dLEx0ixfj;sparse(N-1,N+1);
%                 dLExfix0j,sparse(1,N-1),dLExfixfj];
            dLExx{i,j} = [dLEx0ix0j,sparse(1,N);sparse(N-1,N+1);
                dLExfix0j,sparse(1,N-1),dLExfixfj];
            % Partials of LI wrt xx
            dLIxx{i,j} = [diag(0.5*(tf-t0)*(sigma*w.'.*dgxixj - ...
                sum(reshape(gamma,N,Nx).*daxixj,2)) + sum(reshape(phi,...
                N,Np).*dcxixj,2)),sparse(N,1);sparse(1,N),0];
            if i~=j
%                 dLExx{j,i} = dLExx{i,j}.';
%                 dLIxx{j,i} = dLIxx{i,j}.';
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
            daxiuj = (a(tgr,xgr+hgr.*ex(i,:),u+hu.*eu(jj,:),auxdata)-...
                a(tgr,xgr+hgr.*ex(i,:),u,auxdata) - a(tgr,xgr,...
                u+hu.*eu(jj,:),auxdata)+dyn)./(hgr(:,i).*hu(:,jj));
            if ~isempty(po)
                dgxiuj = (g(tgr,xgr+hgr.*ex(i,:),u+hu.*eu(jj,:),auxdata)-...
                    g(tgr,xgr+hgr.*ex(i,:),u,auxdata)-g(tgr,xgr,u+...
                    hu.*eu(jj,:),auxdata)+po)./(hgr(:,i).*hu(:,jj));
            else
                dgxiuj = zeros(N,1);
            end
            if ~isempty(pc)
                dcxiuj = (c(tgr,xgr+hgr.*ex(i,:),u+hu.*eu(jj,:),auxdata)-...
                    c(tgr,xgr+hgr.*ex(i,:),u,auxdata)-c(tgr,xgr,u+...
                    hu.*eu(jj,:),auxdata)+pc)./(hgr(:,i).*hu(:,jj));
            else
                dcxiuj = zeros(N,1);
            end
            % Partials of LI wrt xu
            dLIxu{jj,i} = [diag(0.5*(tf-t0)*(sigma*w.'.*dgxiuj - ...
                sum(reshape(gamma,N,Nx).*daxiuj,2)) + sum(reshape(phi,...
                N,Np).*dcxiuj,2)),sparse(N,1)];
        end
        
        % Derive wrt Time
        dadx = (a(tgr,xgr+hgr.*ex(i,:),u,auxdata)-...
            a(tgr,xgr-hgr.*ex(i,:),u,auxdata))./(2*hgr(:,i));
        dgdx = (g(tgr,xgr+hgr.*ex(i,:),u,auxdata)-...
            g(tgr,xgr-hgr.*ex(i,:),u,auxdata))./(2*hgr(:,i));
        datxi = (a(tgr+htgr,xgr+hgr.*ex(i,:),u,auxdata)-...
            a(tgr,xgr+hgr.*ex(i,:),u,auxdata) - a(tgr+htgr,xgr,...
            u,auxdata)+dyn)./(hgr(:,i).*htgr);
        if ~isempty(po)
            dgtxi = (g(tgr+htgr,xgr+hgr.*ex(i,:),u,auxdata)-...
                g(tgr,xgr+hgr.*ex(i,:),u,auxdata) - g(tgr+htgr,xgr,...
                u,auxdata)+po)./(hgr(:,i).*htgr);
        else
            dgtxi = zeros(N,1);
        end
        if ~isempty(pc)
            dctxi = (c(tgr+htgr,xgr+hgr.*ex(i,:),u,auxdata)-...
                c(tgr,xgr+hgr.*ex(i,:),u,auxdata) - c(tgr+htgr,xgr,...
                u,auxdata)+pc)./(hgr(:,i).*htgr);
        else
            dctxi = zeros(N,1);
        end
        dLIt0x{1,i} = [(0.5*(sum(reshape(gamma,N,Nx).*dadx,2)-sigma*w.'.*dgdx)+...
            0.5*(tf-t0)*alpha.*(sigma*w.'.*dgtxi-sum(reshape(gamma,N,Nx).*datxi,2))...
            +alpha.*(sum(reshape(phi,N,Np).*dctxi,2))).',0];
        dLItfx{1,i} = [(0.5*(-sum(reshape(gamma,N,Nx).*dadx,2)+sigma*w.'.*dgdx)+...
            0.5*(tf-t0)*beta.*(sigma*w.'.*dgtxi-sum(reshape(gamma,N,Nx).*datxi,2))...
            +beta.*(sum(reshape(phi,N,Np).*dctxi,2))).',0];
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
                dauiuj = (a(tgr,xgr,u+hu.*eu(i,:)+hu.*eu(j,:),auxdata)-...
                    a(tgr,xgr,u+hu.*eu(i,:),auxdata)-...
                    a(tgr,xgr,u+hu.*eu(j,:),auxdata)+dyn)./(hu(:,i).*hu(:,j));
                if ~isempty(po)
                    dguiuj = (g(tgr,xgr,u+hu.*eu(i,:)+hu.*eu(j,:),auxdata)-...
                        g(tgr,xgr,u+hu.*eu(i,:),auxdata)-...
                        g(tgr,xgr,u+hu.*eu(j,:),auxdata)+po)./(hu(:,i).*hu(:,j));
                else
                    dguiuj = zeros(N,1);
                end
                if ~isempty(pc)
                    dcuiuj = (c(tgr,xgr,u+hu.*eu(i,:)+hu.*eu(j,:),auxdata)-...
                        c(tgr,xgr,u+hu.*eu(i,:),auxdata)-...
                        c(tgr,xgr,u+hu.*eu(j,:),auxdata)+pc)./(hu(:,i).*hu(:,j));
                else
                    dcuiuj = zeros(N,1);
                end
            dLIuu{i,j} = diag(0.5*(tf-t0)*(sigma*w.'.*dguiuj - ...
                sum(reshape(gamma,N,Nx).*dauiuj,2)) + sum(reshape(phi,...
                N,Np).*dcuiuj,2));
            if j ~= i
%                 dLIuu{j,i} = dLIuu{i,j}.';
                dLIuu{j,i} = sparse(N,N);
            end
        end
        % Derive wrt Time
        dadu = (a(tgr,xgr,u+hu.*eu(i,:),auxdata)-...
            a(tgr,xgr,u-hu.*eu(i,:),auxdata))./(2*hu(:,i));
        dgdu = (g(tgr,xgr,u+hu.*eu(i,:),auxdata)-...
            g(tgr,xgr,u-hu.*eu(i,:),auxdata))./(2*hu(:,i));
        datui = (a(tgr+htgr,xgr,u+hu.*eu(i,:),auxdata)-...
            a(tgr,xgr,u+hu.*eu(i,:),auxdata) - a(tgr+htgr,xgr,...
            u,auxdata)+dyn)./(hu(:,i).*htgr);
        if ~isempty(po)
            dgtui = (g(tgr+htgr,xgr,u+hu.*eu(i,:),auxdata)-...
                g(tgr,xgr,u+hu.*eu(i,:),auxdata) - g(tgr+htgr,xgr,...
                u,auxdata)+po)./(hu(:,i).*htgr);
        else
            dgtui = zeros(N,1);
        end
        if ~isempty(pc)
            dctui = (c(tgr+htgr,xgr,u+hu.*eu(i,:),auxdata)-...
                c(tgr,xgr,u+hu.*eu(i,:),auxdata) - c(tgr+htgr,xgr,...
                u,auxdata)+pc)./(hu(:,i).*htgr);
        else
            dctui = zeros(N,1);
        end
        dLIt0u{1,i} = (0.5*(sum(reshape(gamma,N,Nx).*dadu,2)-sigma*w.'.*dgdu)+...
            0.5*(tf-t0)*alpha.*(sigma*w.'.*dgtui-sum(reshape(gamma,N,Nx).*datui,2))...
            +alpha.*(sum(reshape(phi,N,Np).*dctui,2))).';
        dLItfu{1,i} = (0.5*(-sum(reshape(gamma,N,Nx).*dadu,2)+sigma*w.'.*dgdu)+...
            0.5*(tf-t0)*beta.*(sigma*w.'.*dgtui-sum(reshape(gamma,N,Nx).*datui,2))...
            +beta.*(sum(reshape(phi,N,Np).*dctui,2))).';
    end
    
    % LI wrt t0
    dadt = (a(tgr+htgr,xgr,u,auxdata)-a(tgr-htgr,xgr,u,auxdata))./(2*htgr);
    dgdt = (g(tgr+htgr,xgr,u,auxdata)-g(tgr-htgr,xgr,u,auxdata))./(2*htgr);
    datt = (a(tgr+2*htgr,xgr,u,auxdata)-2*a(tgr+htgr,xgr,u,auxdata)+...
        dyn)./(htgr.^2);
    if ~isempty(po)
        dgtt = (g(tgr+2*htgr,xgr,u,auxdata)-2*g(tgr+htgr,xgr,u,auxdata)+...
            po)./(htgr.^2);
    else
        dgtt = zeros(N,1);
    end
    if ~isempty(pc)
        dctt = (c(tgr+2*htgr,xgr,u,auxdata)-2*c(tgr+htgr,xgr,u,auxdata)+...
            pc)./(htgr.^2);
    else
        dctt = zeros(N,1);
    end
    
    dLIt0t0 = alpha.'*(sum(reshape(gamma,N,Nx).*dadt,2)-sigma*w.'.*dgdt)+...
        0.5*(tf-t0)*alpha.'*((sigma*w.'.*dgtt-sum(reshape(gamma,N,...
        Nx).*datt,2)).*alpha)+alpha.'*((sum(reshape(phi,N,Np).*dctt,2)).*alpha);
    dLItft0 = 0.5*alpha.'*(-sum(reshape(gamma,N,Nx).*dadt,2)+sigma*w.'.*dgdt)+...
        0.5*beta.'*(sum(reshape(gamma,N,Nx).*dadt,2)-sigma*w.'.*dgdt)+...
        0.5*(tf-t0)*alpha.'*((sigma*w.'.*dgtt-sum(reshape(gamma,N,...
        Nx).*datt,2)).*beta)+alpha.'*((sum(reshape(phi,N,Np).*dctt,2)).*beta);
    dLItftf = beta.'*(-sum(reshape(gamma,N,Nx).*dadt,2)+sigma*w.'.*dgdt)+...
        0.5*(tf-t0)*beta.'*((sigma*w.'.*dgtt-sum(reshape(gamma,N,...
        Nx).*datt,2)).*beta)+beta.'*((sum(reshape(phi,N,Np).*dctt,2)).*beta);
    
    % Partials of LE
%     dLE = [cell2mat(dLExx), dLEux.', cell2mat(dLEt0x).', cell2mat(dLEtfx).';
%         dLEux,           dLEuu,   dLEt0u.',           dLEtfu.'          ;
%         cell2mat(dLEt0x),dLEt0u,  dLEt0t0,            dLEt0tf.'         ;
%         cell2mat(dLEtfx),dLEtfu,  dLEt0tf,            dLEtftf];
        % Partials of LE
    dLE = [cell2mat(dLExx) , sparse(Nx*(N+1),N*Nu), sparse(Nx*(N+1),1),  sparse(Nx*(N+1),1);
            dLEux,            dLEuu,                 sparse(Nu*N,1),     sparse(Nu*N,1);
            cell2mat(dLEt0x), dLEt0u,                dLEt0t0,            sparse(0);
            cell2mat(dLEtfx), dLEtfu,                dLEt0tf,            dLEtftf];
    % Partials of LE
%     dLI = [cell2mat(dLIxx) , cell2mat(dLIxu).', cell2mat(dLIt0x).', cell2mat(dLItfx).';
%         cell2mat(dLIxu) , cell2mat(dLIuu)  , cell2mat(dLIt0u).', cell2mat(dLItfu).';
%         cell2mat(dLIt0x), cell2mat(dLIt0u) , dLIt0t0           , dLItft0.'         ;
%         cell2mat(dLItfx), cell2mat(dLItfu) , dLItft0           , dLItftf           ];
    dLI = [cell2mat(dLIxx) , sparse(Nx*(N+1),N*Nu), sparse(Nx*(N+1),1), sparse(Nx*(N+1),1);
        cell2mat(dLIxu) , cell2mat(dLIuu)  , sparse(Nu*N,1), sparse(Nu*N,1);
        cell2mat(dLIt0x), cell2mat(dLIt0u) , dLIt0t0           , sparse(0)         ;
        cell2mat(dLItfx), cell2mat(dLItfu) , dLItft0           , dLItftf           ];
    rowcolPhase = HessShift+1:HessShift+DecVarSum(iphase);
%     H(rowcolPhase,rowcolPhase) = tril(dLE + dLI);
    H(rowcolPhase,rowcolPhase) = dLE + dLI;
    HessShift = rowcolPhase(end);
end
% H = MAIN.scaling.invV.'*H*MAIN.scaling.invV;
end

function output = LE(t0,tf,x0,xf,auxdata,sigma,nu)
global MAIN
MayerCost = MAIN.funcs.BndObj(t0,tf,x0,xf,auxdata);
BndConstraints = MAIN.funcs.BndCst(t0,tf,x0,xf,auxdata);
output = sigma*MayerCost  + nu*BndConstraints(:);
end