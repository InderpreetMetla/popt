function Jtotal = jacobianstructure()
%*************************************************************%
% This function computes the sparsity structure of the Jacobian 
%
% Outputs: 
%    - Jtotal  : Jacobian Sparsity Pattern
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%*************************************************************%
global MAIN

% Number of Nodes 
Nodes = MAIN.nodes;
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
nlinks = MAIN.nphases-1;
% Allocate Memory for Total Jacobian
Jtotal = sparse(MAIN.nNonLinCsts,sum(MAIN.nPhaseDecVar));
rowshift = 0;
colshiftPhase = 0;
for iphase = 1:MAIN.nphases
    N   = Nodes(iphase);
    Nx  = nStates(iphase);
    Nu  = nControls(iphase);
    Np  = nPathCsts(iphase);
    Nb  = nBndCsts(iphase);
    DefectsparsityPerState = MAIN.LGColloc.phase(iphase).DM1;
    row = rowshift+1:rowshift+MAIN.nPhaseCstsLessBnd(iphase);
    colshift = colshiftPhase;
    col = colshift+1:colshift+Nx*(N+1); 
    Jtotal(row,col) = [kron(eye(Nx),DefectsparsityPerState)+...
        repmat(sparse([diag(ones(N,1)),zeros(N,1)]),Nx,Nx);...
        repmat(sparse([diag(ones(N,1)),zeros(N,1)]),Np,Nx)];
    colshift = col(end);
    col = colshift+1:colshift+N*Nu;
    Jtotal(row,col) = repmat(diag(ones(N,1)),Nx+Np,Nu);
    colshift = col(end);
    col = colshift+1:colshift+2;
    Jtotal(row,col) = [ones(N*Nx,2);...
        [ones(N*Np,1),ones(N*Np,1)]];
    rowshift = row(end);
    colshift = colshiftPhase;
    row = rowshift+1:rowshift+Nb;
    col = colshift+[(1:Nx)*(N+1)-N,(1:Nx)*(N+1),N*(Nx+Nu)+Nx+[1:2]];
    % Boundary wrt Boundary States
    Jtotal(row,col) = 1;
    rowshift = rowshift+Nb;
    colshiftPhase = col(end);
end
if nlinks > 0
    Jtotal(MAIN.LinkJac.rowIdx,:) = spones(MAIN.LinkJac.Jlinear);
end
Jtotal = spones(Jtotal);
end