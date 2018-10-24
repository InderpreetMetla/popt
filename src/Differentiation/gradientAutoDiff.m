function g = gradientAutoDiff(x)
%*********************************************************************%
% This function uses the auto differentiation package MatlabAutoDiff
% by martinReasearch to compute the gradient of the objective function
%
% Inputs: 
%    - x   : Vector of decision variables
%
% Outputs: 
%    - g   : gradient of objFun.m
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%*********************************************************************%
g = AutoDiffJacobianAutoDiff(@objFun,x);
end