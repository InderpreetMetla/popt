function funcs = ipoptfuncSetupFD()
%******************************************************************%
% This function sets up the funcs struct for ipopt if the user has 
% request First Order Finite Difference Differentiation
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
funcs.gradient          = @gradientFD;
funcs.jacobian          = @jacobianFD;
funcs.jacobianstructure = @jacobianstructure;
end