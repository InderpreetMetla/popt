function funcs = DerivMethod(problem)
%******************************************************************%
% This function assigns the specific derivation method requested
% by the user for ipopt. If no method specified then FD is used.
%
% Inputs:
% - problem :	User supplied problem struct + new fields
%               from previous functions
%
% Outputs:
%  - funcs:     struct with following fields
%               - funcs.objective
%               - funcs.constraints
%               - funcs.gradient
%               - funcs.jacobian
%               - funcs.jacobianstructure
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%******************************************************************%

if ~isfield(problem.derivatives,'method') || ...
        isequal(problem.derivatives.method,'fd') || ...
        isequal(problem.derivatives.method,'FD') || ...
        isequal(problem.derivatives.method,'forward-difference')
    
    if isfield(problem.derivatives,'order') && ...
            (isequal(problem.derivatives.order,2) || ...
            isequal(problem.derivatives.order,'2'))
        funcs = ipoptfuncSetupFDNewton(problem.autoscale);
        HessFlag = 1;
    else
        funcs = ipoptfuncSetupFD;
        HessFlag = 0;
    end
    
    if problem.grid.iter == 0 && HessFlag == 1
        disp('Calculating first and second derivatives using Forward Differences...')
    elseif problem.grid.iter == 0
        disp('Calculating first derivatives using Forward Differences...')
    end
    
elseif isequal(problem.derivatives.method,'CS') || ...
        isequal(problem.derivatives.method,'cs')
    
    if isfield(problem.derivatives,'order') && ...
            (isequal(problem.derivatives.order,2) || ...
            isequal(problem.derivatives.order,'2'))
        funcs = ipoptfuncSetupCSNewton(problem.autoscale);
        HessFlag = 1;
    else
        funcs = ipoptfuncSetupCS;
        HessFlag = 0;
    end
    
    if problem.grid.iter == 0 && HessFlag == 1
        disp('Calculating first derivatives using Complex Step Differentiation and second derivatives using Forward Differences...')
    elseif problem.grid.iter == 0
        disp('Calculating first derivatives using Complex Step Differentiation...')
    end
    
elseif isequal(problem.derivatives.method,'AD') || ...
        isequal(problem.derivatives.method,'ad')
    
    if isfield(problem.derivatives,'order') && ...
            (isequal(problem.derivatives.order,2) || ...
            isequal(problem.derivatives.order,'2'))
        funcs = ipoptfuncSetupADNewton(problem.autoscale);
        HessFlag = 1;
    else
        funcs = ipoptfuncSetupAD;
        HessFlag = 0;
    end
    
    if problem.grid.iter == 0 && HessFlag == 1
        disp('Calculating first derivatives using Automatic Differentiation and second derivatives using Forward Differences...')
    elseif problem.grid.iter == 0
        disp('Calculating first derivatives using Automatic Differentiation...')
    end
else
    error(strcat('Error: derivative.method must be ''FD'',',...
        ' ''CS'' or ''AD'''));
end


end