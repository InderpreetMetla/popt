%******************************************************************%
% Infinite Horizon LQR from:
% Kirk, D. E. (2012). “Optimal control theory: an introduction”. 
% Courier Corporation, pp. 216-217. 
%
% minimise: 
%       \int_0^\infty (x_1^2+ 0.5*x_2^2+0.25*u^2) dt
% subject to:
%       x_1dot = x_2
%       x_2dot = 2*x_1 - x_2 + u
% intitial condition:
%       x_1(0) = -4 and x_2(0) = 4
%******************************************************************%

clear;clc

t0 = -1; 
tf = 1;
x1_0 = -4; x2_0 = 4;
x1min = -10;
x1max = -x1min;
x2min = -10;
x2max = -x2min;
xf = 1;

umin = -100;
umax = -umin;
%******************************************************************%
%  Setup for Problem Bounds
%******************************************************************%
bounds.phase.initialtime.lb     = t0; 
bounds.phase.initialtime.ub     = t0;
bounds.phase.finaltime.lb       = tf; 
bounds.phase.finaltime.ub       = tf;

bounds.phase.initialstate.lb    = [x1_0, x2_0]; 
bounds.phase.initialstate.ub    = [x1_0, x2_0]; 
bounds.phase.state.lb           = [x1min, x2min];
bounds.phase.state.ub           = [x1max, x2max];
bounds.phase.finalstate.lb      = [x1min, x2min];
bounds.phase.finalstate.ub      = [x1max, x2max];

bounds.phase.control.lb         = umin; 
bounds.phase.control.ub         = umax;

%******************************************************************%
% Provide Guess of Solution  
%******************************************************************% 
guess.phase.time    = [t0; tf]; 
guess.phase.state   = [[x1_0;x1max],[x2_0;x2max]];
guess.phase.control = [umin;umax];

%******************************************************************%
% Setup Grid
%******************************************************************% 
grid.phase.nodes.initialgrid = [5,5,5,5,3,2,2];
grid.phase.nodes.lb = 4;
grid.phase.nodes.ub = 12;
grid.max_refine     = 50;
grid.tol            = 1e-6;

%******************************************************************%
% Assemble Information into Problem Structure  
%******************************************************************%
problem.name                = 'lqrIH-Problem';
problem.funcs.Dynamics      = @lqrIHDynamics; 
problem.funcs.PathObj       = @lqrIHPathObj;
problem.funcs.BndObj        = @lqrIHBndObj; 
problem.funcs.PathCst       = @lqrIHPathCst;
problem.funcs.BndCst        = @lqrIHBndCst;
problem.bounds              = bounds;
problem.guess               = guess;
problem.derivatives.method	= 'cs';
problem.derivatives.order   = 2;
problem.grid                = grid;
problem.autoscale           = 'off';
problem.options.ipopt.linear_solver = 'ma57';

%******************************************************************%
% Solve Problem
%******************************************************************%
[solution, plotdata] = popt(problem);