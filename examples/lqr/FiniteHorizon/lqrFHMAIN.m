%******************************************************************%
% Finite Horizon LQR
%******************************************************************%

clear;clc

t0 = 0; 
tf = 5;
x0 = 2;
xmin = -10;
xmax = -xmin;
xf = 1;

umin = -10;
umax = -umin;
%******************************************************************%
%  Setup for Problem Bounds
%******************************************************************%
bounds.phase.initialtime.lb     = t0; 
bounds.phase.initialtime.ub     = t0;
bounds.phase.finaltime.lb       = tf; 
bounds.phase.finaltime.ub       = tf;

bounds.phase.initialstate.lb    = x0; 
bounds.phase.initialstate.ub    = x0; 
bounds.phase.state.lb           = xmin; 
bounds.phase.state.ub           = xmax; 
bounds.phase.finalstate.lb      = xf; 
bounds.phase.finalstate.ub      = xf; 

bounds.phase.control.lb         = umin; 
bounds.phase.control.ub         = umax;

%******************************************************************%
% Provide Guess of Solution  
%******************************************************************% 
guess.phase.time    = [t0; tf]; 
guess.phase.state   = [x0;xf];
guess.phase.control = [umin;umax];

%******************************************************************%
% Setup Grid
%******************************************************************% 
grid.phase.nodes.initialgrid = 3;
grid.phase.nodes.lb = 3;
grid.phase.nodes.ub = 10;
grid.tol            = 1e-6;

%******************************************************************%
% Assemble Information into Problem Structure  
%******************************************************************%
problem.name                = 'lqrFH-Problem';
problem.funcs.Dynamics      = @lqrFHDynamics; 
problem.funcs.PathObj       = @lqrFHPathObj;
problem.funcs.BndObj        = @lqrFHBndObj; 
problem.funcs.PathCst       = @lqrFHPathCst;
problem.funcs.BndCst        = @lqrFHBndCst;
% problem.auxdata             = auxdata;
problem.bounds              = bounds;
problem.guess               = guess;
problem.derivatives.method	= 'fd';
problem.derivatives.order    = 2;
problem.grid                = grid;
problem.options.ipopt.linear_solver = 'ma57';

%******************************************************************%
% Solve Problem
%******************************************************************%
tic
[solution, plotdata] = popt(problem);
toc

%******************************************************************%
% Plotting
%******************************************************************%
t = plotdata.phase.time;
x = plotdata.phase.state;
u = plotdata.phase.control;

plot(t,x(:,1),'-o')
hold on
plot(t,u,'-o')
xlabel('Time')
ylabel('Solution')