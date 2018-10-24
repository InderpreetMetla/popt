function funcs = ipoptfuncSetupCSNewton(autoscale)
%******************************************************************%
% This function sets up the funcs struct for ipopt if the user has 
% request Complex Step differentiation. Second derivatives are 
% computed using Finite Differences. 
%
% Outputs:
%  - funcs:     struct with following fields
%               - funcs.objective
%               - funcs.constraints
%               - funcs.gradient
%               - funcs.jacobian
%               - funcs.jacobianstructure
%               - funcs.hessian
%               - funcs.hessianstructure
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%******************************************************************%
funcs.objective         = @objFun;
funcs.constraints       = @cstFun;
funcs.gradient          = @gradientCS;
funcs.jacobian          = @jacobianCS;
funcs.jacobianstructure = @jacobianstructure;
if isequal(autoscale,'on') || isequal(autoscale,1)
    funcs.hessian       = @hessianFDScaled;
else
    funcs.hessian       = @hessianFD;
end
funcs.hessianstructure  = @hessianstructure;
end