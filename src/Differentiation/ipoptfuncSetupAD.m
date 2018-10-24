function funcs = ipoptfuncSetupAD()
%******************************************************************%
% This function sets up the funcs struct for ipopt if the user has 
% request First Order Automatic Differentiation using MatlabAutoDiff
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
funcs.gradient          = @gradientAutoDiff;
funcs.jacobian          = @jacobianAutoDiff;
funcs.jacobianstructure = @jacobianstructure;
end