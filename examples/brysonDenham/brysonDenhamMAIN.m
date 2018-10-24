%******************************************************************%
% Bryson Denham Problem 
%******************************************************************%

clear;%clc

L       = 1/9;
t0      = 0;     
tf      = 1;
x0  	= 0;
xf    	= 0;
xmin  	= 0;
xmax  	= L;
v0      = 1;
vf      = -v0;
vmin	= -10;
vmax	= 10;
umin	= -10;
umax	= 5;

%******************************************************************%
%  Setup for Problem Bounds  
%******************************************************************%
bounds.phase.initialtime.lb     = t0; 
bounds.phase.initialtime.ub     = t0;
bounds.phase.finaltime.lb       = tf; 
bounds.phase.finaltime.ub       = tf;

bounds.phase.initialstate.lb    = [x0,v0]; 
bounds.phase.initialstate.ub    = [x0,v0]; 
bounds.phase.state.lb           = [xmin,vmin]; 
bounds.phase.state.ub           = [vmax,vmax]; 
bounds.phase.finalstate.lb      = [xf,vf]; 
bounds.phase.finalstate.ub      = [xf,vf]; 

bounds.phase.control.lb         = umin; 
bounds.phase.control.ub         = umax;

bounds.phase.path.lb            = xmin; 
bounds.phase.path.ub            = xmax;

%******************************************************************%
% Provide Guess of Solution  
%******************************************************************%
guess.phase.time    = [t0; tf]; 
guess.phase.state   = [[x0; xf],[v0; vf]];
guess.phase.control = [umin; 0];
grid.tol = 1e-10;
grid.phase.nodes.lb = 4;
grid.phase.nodes.ub = 10;
grid.phase.nodes.initialgrid = 8*ones(1,4);

%******************************************************************%
% Assemble Information into Problem Structure  
%******************************************************************%
problem.name                = 'brysonDenham';
problem.funcs.Dynamics      = @brysonDenhamDynamics;
problem.funcs.PathObj       = @brysonDenhamPathObj;
problem.funcs.PathCst       = @brysonDenhamPathCst;
problem.funcs.BndObj        = @brysonDenhamBndObj; 
problem.funcs.BndCst        = @brysonDenhamBndCst;
problem.bounds              = bounds;
problem.guess               = guess;
problem.derivatives.method	= 'fd';
problem.derivatives.first.stepsize.jacobian	= 1.9e-5;
problem.derivatives.order    = 1;
problem.autoscale           = 'off';
problem.grid                = grid;

%******************************************************************%
% Solve Problem
%******************************************************************%
tic
[solution, plotdata]   = popt(problem);
toc

%******************************************************************%
% Plotting
%******************************************************************%
t = solution.phase.time;
x = solution.phase.state;
u = solution.phase.control;
[x_exact,v_exact,u_exact,t_exact] = brysonDenhamExactSol(L,t);

subplot(1,3,1)
plot(t,x(:,1),'ko-','MarkerSize',10)
hold on
plot(t_exact,x_exact,'k-','LineWidth',2)
xlabel('Time ($s$)','Interpreter','LaTeX','FontSize',15)
ylabel('Position ($m$)','Interpreter','LaTeX','FontSize',15)
legend({'POPT', 'Exact'},'Interpreter','LaTeX','FontSize',16)
subplot(1,3,2)
plot(t,x(:,2),'ko-','MarkerSize',10)
hold on
plot(t_exact,v_exact,'k-','LineWidth',2)
xlabel('Time ($s$)','Interpreter','LaTeX','FontSize',15)
ylabel('Velocity ($m/s$)','Interpreter','LaTeX','FontSize',15)
subplot(1,3,3)
plot(t,u,'ko-','MarkerSize',10)
hold on
plot(t_exact,u_exact,'k-','LineWidth',2)
xlabel('Time ($s$)','Interpreter','LaTeX','FontSize',15)
ylabel('Control ($m/s^2$)','Interpreter','LaTeX','FontSize',15)
% max_error.x = max(abs(x(:,1)-x_exact.'));
% max_error.v = max(abs(x(:,2)-v_exact.'));
% max_error.u = max(abs(u-u_exact.'));
% 
% max_error
% figure(2)
% for i=1:length(solution.grid_analysis.iteration)
%     plot(solution.grid_analysis.iteration(i).phase.node_locations.',i,'ko')
%     hold on
% end