%******************************************************************%
% Hyper-Sensitive Problem
%******************************************************************%
% This example is taken from the following reference:              
% Rao, A. V., and Mease, K. D., "Eigenvector Approximate Dichotomic
% Basis Methods for Solving Hyper-Sensitive Optimal Control 
% Problems," Optimal Control Applications and Methods, Vol. 21, 
% No. 1., January-February 2000, pp. 1-17.          
%******************************************************************%

clear;clc
t0 = 0;
tf = 10000;
x0 = 1;
xf = 1.5; 
xMin = -500;
xMax = +500;
uMin = -500;
uMax = +500;

%******************************************************************%
%  Setup for Problem Bounds
%******************************************************************%
bounds.phase.initialtime.lb = t0;
bounds.phase.initialtime.ub = t0;
bounds.phase.finaltime.lb = tf;
bounds.phase.finaltime.ub = tf;
bounds.phase.initialstate.lb = x0; 
bounds.phase.initialstate.ub = x0;
bounds.phase.state.lb = xMin;
bounds.phase.state.ub = xMax;
bounds.phase.finalstate.lb = xf; 
bounds.phase.finalstate.ub = xf;
bounds.phase.control.lb = uMin; 
bounds.phase.control.ub = uMax;

%******************************************************************%
% Provide Guess of Solution  
%******************************************************************% 
guess.phase.time     = [t0; tf]; 
guess.phase.state    = [x0; xf];
guess.phase.control  = [0; 0];

%******************************************************************%
% Setup Grid
%******************************************************************% 
grid.phase.nodes.initialgrid = 8;
grid.phase.nodes.lb = 2;
grid.phase.nodes.ub = 14;
grid.tol            = 1e-7;
grid.max_refine     = 10;

%******************************************************************%
% Setup Options
%******************************************************************% 
options.ipopt.linear_solver = 'ma57';

%******************************************************************%
% Assemble Information into Problem Structure  
%******************************************************************%
problem.name                = 'hypersensitive-Problem';
problem.funcs.Dynamics      = @hypersensitiveDynamics; 
problem.funcs.PathObj       = @hypersensitivePathObj;
problem.funcs.BndObj        = @hypersensitiveBndObj; 
problem.funcs.PathCst       = @hypersensitivePathCst;
problem.funcs.BndCst        = @hypersensitiveBndCst;
problem.bounds              = bounds;
problem.guess               = guess;
problem.derivatives.method	= 'fd';
problem.derivatives.order    = 2;
problem.autoscale           = 'on';
problem.grid                = grid;
problem.options             = options;

%******************************************************************%
% Solve Problem
%******************************************************************%
[solution, plotdata] = popt(problem);
Total_Time = solution.cpu_time.total;
%******************************************************************%
% Plotting
%******************************************************************%
t = solution.phase.time;
x = solution.phase.state;
u = solution.phase.control;
hold on
subplot(2,1,1)
plot(t,x(:,1),'-o')
xlabel('Time')
ylabel('Horizontal Position')
hold on
subplot(2,1,2)
plot(t,u,'-o')
xlabel('Time')
ylabel('Control')