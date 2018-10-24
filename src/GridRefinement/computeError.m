function RelativeError = computeError(solution,auxdata,TauOld,...
        NodesPerInt,BreakPts,nInts,iphase,DynFunc,stateScales)
%**********************************************************************%
% This function computes the relative error in the solution.
% It does this by creating a solution on N+1 points and then integrating
% the dynamics using the Radau Integration Matrix. This is compared
% to an interpolated solution fvrom ipopt. 
%
% Inputs:
%
% - solution    : The solution from ipopt
% - auxdata     : Auxilary data for the problem
% - TauOld      : The collocation points on the grid that the current
%                 solution was obtained
% - NodesPerInt : No. of nodes per interval on current grid
% - BreakPts    : Location of break points on current grid
% - nInts       : No. of intervals on current grid
% - iphase      : Phase number
% - DynFunc     : Function handle for the Dynamics function
% - stateScales : The state scales to scale the dynamics function
%
% Outputs:
%
% - RelativeError : The relative error for this grid
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%**********************************************************************%
auxdata.iphase = iphase;
Idx = [1 cumsum(NodesPerInt)+ones(1,nInts)];
Tau = []; X = Tau; LagrangeControl = Tau;

for iInt = 1:nInts
    N       = NodesPerInt(iInt);
    t0Int   = TauOld(Idx(iInt));
    tfInt   = TauOld(Idx(iInt+1));
    Xold    = solution.phase(iphase).state(Idx(iInt):Idx(iInt+1),:);
    TauNew  = legsrd(N+1);
    TauNew	= [0.5*(tfInt-t0Int).*(TauNew+1)+t0Int;BreakPts(iInt+1)];
    XNew    = barylag([TauOld(Idx(iInt):Idx(iInt+1)),...
        Xold],TauNew);
    if ~isempty(solution.phase(iphase).control)
        Uold = solution.phase(iphase).control(Idx(iInt):...
            Idx(iInt+1),:);
        UNew = barylag([TauOld(Idx(iInt):...
            Idx(iInt+1)-1),Uold(1:end-1,:)],TauNew(1:end-1));
    else
        UNew = [];
    end
    Tau = [Tau;TauNew(1:end-1)];
    X = [X;XNew(1:end-1,:)];
    LagrangeControl = [LagrangeControl;UNew];
end

t0 = solution.phase(iphase).time(1);
tf = solution.phase(iphase).time(end);
LagrangeTime = 0.5 * (tf - t0) .* ([Tau;1] + 1) + t0;
LagrangeState = [X;solution.phase(iphase).state(end,:)];

F = DynFunc(LagrangeTime(1:end-1),LagrangeState(1:end-1,:),...
    LagrangeControl,auxdata);
Gauss = GaussRadauMethod(NodesPerInt+1,BreakPts,size(NodesPerInt+1,2));
Yhat0 = Gauss.MM*LagrangeState;
Yhat = Yhat0 + 0.5*(tf-t0)*Gauss.IM*F;
AbsError = abs([LagrangeState(1,:);Yhat]-LagrangeState);
AbsError = AbsError./stateScales{iphase};
% RelativeError = AbsError./(1+max(abs(LagrangeState)));
RelativeError = AbsError./(1+max(max(abs(LagrangeState)),max(abs(F))));
end