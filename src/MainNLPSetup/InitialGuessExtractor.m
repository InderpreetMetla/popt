function problem = InitialGuessExtractor(problem)
%******************************************************************%
% This function extracts the initial guess that is provided by 
% the user.
% 
% The guess is transformed to exist over a closed time period 
% [-1,1] and the Gauss Radau Pseudospectral method (GRPM) is used   
% to determine the collocation points, i.e. where the dynamics    
% and constraints will be satisfied. 
% The initial guess is interpolated at every collocation point.
% 
% Inputs:
% - problem :	User supplied problem problem struct + new fields  
%               from previous functions
%
% Outputs:
%  - problem struct with:
%       - .DecVar0 field :  Initial Guess for ipopt
%       - .LGColloc      :  struct of GRPM differentiation matrix,
%                           collocation pts and quadrature weights
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%******************************************************************%

guess       = problem.guess.phase;
nodes       = problem.nodes;
nphases     = problem.nphases;
sizes       = problem.sizes;
nStates     = sizes(1,:);
nControls   = sizes(2,:);
% Assign an empty cell for initial guess vector
DecVar0     = cell(nphases,1);

for iphase=1:nphases
    
    % Extract Guess
    tGuess      = guess(iphase).time;
    stateGuess  = guess(iphase).state;
    
    if isfield(guess(iphase),'control')
        controlGuess = guess(iphase).control;
    else
        controlGuess = [];
    end
    
    % Interpolate guess
    NodesPerInt = problem.grid.phase(iphase).nodes.PerInterval;
    
    BreakPts    = problem.grid.phase(iphase).BreakPts; 
    NumInts     = size(NodesPerInt,2);
    LGR         = GaussRadauMethod(NodesPerInt,...
        BreakPts,NumInts);
    Tau_All_Points  = [LGR.Points; 1];
    
    t0Guess         = tGuess(1);
    tfGuess         = tGuess(end);
    
    if isequal(t0Guess,tfGuess)
        error(strcat('Error: Difference between the guess for the',...
            ' Initial and Final Time in Phase %i must',...
            ' be >= 1e-6.'),iphase);
    end

    timeInterpolation = affineTimeTransform(tGuess,Tau_All_Points);
    
    if nStates(iphase)>0
        stateInterpolation = GRPMInterp(tGuess, ...
            stateGuess,timeInterpolation,'state');
    else
        stateInterpolation = [];
    end
    
    if nControls(iphase)>0
        controlInterpolation = GRPMInterp(tGuess, ...
            controlGuess,timeInterpolation,'control');
    else
        controlInterpolation = [];
    end
    
    DecVar0{iphase,1} = [stateInterpolation(:); ...
        controlInterpolation(:); t0Guess; tfGuess];
    
    % Store Phase Specific Legendre Gauss Collocation Parameters
    LGRColloc.phase(iphase).DM       = LGR.DM;
    LGRColloc.phase(iphase).DM1      = LGR.DM1;
    LGRColloc.phase(iphase).Points   = LGR.Points;
    LGRColloc.phase(iphase).Weights  = LGR.Weights;
end

% Output initial guess vector and the Gauss Legendre Pseudospectral
% Collocation struct for and differentiation
problem.DecVar0	 = vertcat(DecVar0{:});
problem.LGColloc = LGRColloc;

end

function Time_Interp = affineTimeTransform(tGuess,tau)
%******************************************************************%
% The function conducts an affine transformation mapping 
% the collocation points from the closed and bounded 
% [-1,1] space to the [t0,tf] space.
%******************************************************************%
t0Guess 	= tGuess(1);
tfGuess     = tGuess(end);
Time_Interp	= 0.5 * (tfGuess - t0Guess) .* (tau + 1) + t0Guess;
end

function Yinterp = GRPMInterp(x,y,Xinterp,Case)
%******************************************************************%
% This function interpolates to find Yinterp = interp1(x,y,Xinterp)
% ie the values of Yinterp at Xinterp of the function y = F(x).
% Xinterp will always be the GPM points plus [-1,1] which have been
% affine transformed to the [t0,tf] domain. 
%
% Case is a switch between state and control
%******************************************************************%
switch Case
    
    case{'state'}
        Yinterp = interp1(x,y,Xinterp,'pchip');
        
    case{'control'}
        Yinterp = interp1(x,y,Xinterp(1:end-1),'pchip');
        
end
end
