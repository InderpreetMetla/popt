function funcs = ipoptfuncSetupCS()
%******************************************************************%
% This function sets up the funcs struct for ipopt if the user has 
% request Complex Step differentiation
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
funcs.objective         = @objFun;
funcs.constraints       = @cstFun;
funcs.gradient          = @gradientCS;
funcs.jacobian          = @jacobianCS;
funcs.jacobianstructure = @jacobianstructure;
end