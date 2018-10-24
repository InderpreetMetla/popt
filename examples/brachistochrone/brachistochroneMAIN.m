%******************************************************************%
% Classical Brachistochrone Problem 
%******************************************************************%
% The problem solved here is given as follows: 
%   Minimize t_f                               
% subject to the dynamic constraints                 
%    dx/dt = v*sin(u)                             
%    dy/dt = v*cos(u)                             
%    dv/dt = g*cos(u)                               
% and the boundary conditions                       
%    x(0) = 0, y(0) = 0, v(0) = 0                  
%    x(t_f) = 2, y(t_f) = 2, v(t_f) = FREE          
%******************************************************************%

clear;clc

auxdata.g = 10;
t0 = 0; 
tfmin = 0; 
tfmax = 10;
x0 = 0; 
y0 = 0; 
v0 = 0;
xf = 2; 
yf = 2;
xmin = 0; 
xmax = 10;
ymin = 0; 
ymax = 10;
vmin = -50; 
vmax = 50;
umin = -pi/2;
umax = pi/2;

%******************************************************************%
%  Setup for Problem Bounds
%******************************************************************%
bounds.phase.initialtime.lb     = t0; 
bounds.phase.initialtime.ub     = t0;
bounds.phase.finaltime.lb       = tfmin; 
bounds.phase.finaltime.ub       = tfmax;

bounds.phase.initialstate.lb    = [x0,y0,v0]; 
bounds.phase.initialstate.ub    = [x0,y0,v0]; 
bounds.phase.state.lb           = [xmin,ymin,vmin]; 
bounds.phase.state.ub           = [xmax,ymax,vmax]; 
bounds.phase.finalstate.lb      = [xf,yf,vmin]; 
bounds.phase.finalstate.ub      = [xf,yf,vmax]; 

bounds.phase.control.lb         = umin; 
bounds.phase.control.ub         = umax;

%******************************************************************%
% Provide Guess of Solution  
%******************************************************************% 
guess.phase.time    = [t0; tfmax]; 
guess.phase.state   = [[x0; xf],[y0; yf],[v0; v0]];
guess.phase.control = [0; pi/2];

%******************************************************************%
% Setup Grid
%******************************************************************% 
grid.phase.nodes.initialgrid = 8;
grid.phase.nodes.lb = 3;
grid.phase.nodes.ub = 10;
grid.tol            = 1e-7;

%******************************************************************%
% Assemble Information into Problem Structure  
%******************************************************************%
problem.name                = 'Brachistochrone-Problem';
problem.funcs.Dynamics      = @brachistochroneDynamics; 
problem.funcs.PathObj       = @brachistochronePathObj;
problem.funcs.BndObj        = @brachistochroneBndObj; 
problem.funcs.PathCst       = @brachistochronePathCst;
problem.funcs.BndCst        = @brachistochroneBndCst;
problem.auxdata             = auxdata;
problem.bounds              = bounds;
problem.guess               = guess;
problem.derivatives.method	= 'fd';
problem.derivatives.order    = 1;
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

subplot(2,2,1)
plot(t,x(:,1),'o')
xlabel('Time')
ylabel('Horizontal Position')

subplot(2,2,2)
plot(t,x(:,2),'o')
xlabel('Time')
ylabel('Vertical Position')

subplot(2,2,3)
plot(t,x(:,3),'o')
xlabel('Time')
ylabel('Velocity')

subplot(2,2,4)
plot(t,u,'o')
xlabel('Time')
ylabel('Control')