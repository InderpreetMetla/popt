function Cost = objFun(Z)
%******************************************************************%
% This function computes the Cost and is called by IPOPT 
% Output is a scalar.
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%******************************************************************%
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
auxdata = MAIN.auxdata;
%*************************************************************%
% Initialise the input for the cost calculation
Cost = 0;
% Unscale
Z = (Z-MAIN.scaling.decvar_shifts)./MAIN.scaling.decvar_scales;

for iphase = 1:nphases
    
    t0      = Z(Idx_xut.phase(iphase).timeIdx(1));
    tf      = Z(Idx_xut.phase(iphase).timeIdx(2));
    Tau     = [LGColloc.phase(iphase).Points;1];
    Time	= 0.5 * (tf - t0) .* (Tau + 1) + t0;
    TimeGR	= Time(1:end-1);
    
    StateVector	= Z(Idx_xut.phase(iphase).stateIdx);
    StateMatrix	= reshape(StateVector,nodes(iphase)+1,nStates(iphase));
    StateGRMatrix = StateMatrix(1:end-1,:);
    x0 = StateMatrix(1,:);
    xf = StateMatrix(end,:);
    
    ControlVector = Z(Idx_xut.phase(iphase).controlIdx);
    
    if nControls(iphase)>0
        ControlMatrix = reshape(ControlVector,...
            nodes(iphase),nControls(iphase));
    else
        ControlMatrix = [];
    end
     
    
    %**************************************************************%
    % Cost Calculation
    %**************************************************************%
    auxdata.iphase = iphase;
    % Compute Objective Cost
    Lagrangian  = MAIN.funcs.PathObj(TimeGR,StateGRMatrix,ControlMatrix,...
        auxdata);
    Mayer       = MAIN.funcs.BndObj(t0,tf,x0,xf,auxdata);
    LGWeights	= LGColloc.phase(iphase).Weights;
    Quadrature	= 0.5 * (tf - t0) * LGWeights * Lagrangian;
    Cost        = Cost + Mayer + Quadrature;
end
    
end