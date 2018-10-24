%******************************************************************%
% Lee-Ramirez bioreactor
%******************************************************************%
        
%******************************************************************%

clear;clc

t0 = 0;
tf = 10;
x10min = 1; x1min = 0; x1fmin = 0;
x20min = 0.1; x2min = 0; x2fmin = 0;
x30min = 40; x3min = 0; x3fmin = 0;
x40min = 0; x4min = 0; x4fmin = 0;
x50min = 0; x5min = 0; x5fmin = 0;
x60min = 1; x6min = 0; x6fmin = 0;
x70min = 0; x7min = 0; x7fmin = 0;
x80min = 0; x8min = 0; x8fmin = 0;
x90min = 0; x9min = 0; x9fmin = 0;

x10max = 1; x1max = 8; x1fmax = 8;
x20max = 0.1; x2max = 8; x2fmax = 8;
x30max = 40; x3max = 45; x3fmax = 45;
x40max = 0; x4max = 8; x4fmax = 8;
x50max = 0; x5max = 8; x5fmax = 8;
x60max = 1; x6max = 8; x6fmax = 8;
x70max = 0; x7max = 8; x7fmax = 8;
x80max = 1; x8max = 1; x8fmax = 1;
x90max = 1; x9max = 1; x9fmax = 1;

u1min = -10;
u1max = 10;
u2min = -10;
u2max = 10;

%******************************************************************%
%  Setup for Problem Bounds
%******************************************************************%
bounds.phase.initialtime.lb     = t0; 
bounds.phase.initialtime.ub     = t0;
bounds.phase.finaltime.lb       = tf; 
bounds.phase.finaltime.ub       = tf; 

bounds.phase.initialstate.lb    = [x10min,x20min,x30min,x40min,x50min,x60min,x70min,x80min,x90min]; 
bounds.phase.initialstate.ub    = [x10max,x20max,x30max,x40max,x50max,x60max,x70max,x80max,x90max]; 
bounds.phase.state.lb           = [x1min,x2min,x3min,x4min,x5min,x6min,x7min,x8min,x9min];
bounds.phase.state.ub           = [x1max,x2max,x3max,x4max,x5max,x6max,x7max,x8max,x9max]; 
bounds.phase.finalstate.lb      = [x1fmin,x2fmin,x3fmin,x4fmin,x5fmin,x6fmin,x7fmin,x8fmin,x9fmin];
bounds.phase.finalstate.ub      = [x1fmax,x2fmax,x3fmax,x4fmax,x5fmax,x6fmax,x7fmax,x8fmax,x9fmax]; 

bounds.phase.control.lb         = [u1min,u2min]; 
bounds.phase.control.ub         = [u1max,u2max];

%******************************************************************%
% Provide Guess of Solution  
%******************************************************************% 
guess.phase.time    = [t0; tf]; 
guess.phase.state   = [[1;4],[0.1;7],[40;40],[0;1],[0;1],[1;1],[0;0],[0;1],[0;1]];
guess.phase.control = [[0;0],[0;0]];

%******************************************************************%
% Setup Grid
%******************************************************************% 
grid.phase.nodes.initialgrid = [20, 35, 50];
grid.phase.nodes.lb = 4;
grid.phase.nodes.ub = 12;
grid.tol            = 1e-3;
grid.max_iter       = 15;

%******************************************************************%
% Assemble Information into Problem Structure  
%******************************************************************%
problem.name                = 'leeBioReactor-Problem';
problem.funcs.Dynamics      = @leeBioReactorDynamics; 
problem.funcs.PathObj       = @leeBioReactorPathObj;
problem.funcs.BndObj        = @leeBioReactorBndObj; 
problem.funcs.PathCst       = @leeBioReactorPathCst;
problem.funcs.BndCst        = @leeBioReactorBndCst;
problem.bounds              = bounds;
problem.guess               = guess;
problem.derivatives.method	= 'ad';
problem.derivatives.order    = 2;
problem.grid                = grid;
problem.autoscale           = 'off';
problem.options.ipopt.linear_solver = 'mumps';

%******************************************************************%
% Solve Problem
%******************************************************************%
tic
[solution, plotdata] = popt(problem);
toc
%%
%******************************************************************%
% Plotting
%******************************************************************%
t = solution.phase.time;
x = solution.phase.state;
u = solution.phase.control;
figure(1)
plot(t,x(:,1),'r-','LineWidth',2)
hold on
plot(t,x(:,2),'g-','LineWidth',2)
plot(t,x(:,3)/10,'b-','LineWidth',2)
plot(t,x(:,4),'m-','LineWidth',2)
plot(t,x(:,5),'c-','LineWidth',2)
plot(t,x(:,6),'-','color',[.48 0.1 0.1],'LineWidth',2)
plot(t,x(:,7),'k-','LineWidth',2)
legend({'$x_1$','$x_2$','$x_3/10$','$x_4$','$x_5$','$x_6$','$x_7$'},...
    'Interpreter','LaTeX','FontSize',18,'Location','NorthWest')
xlabel('Time (s)','Interpreter','LaTeX','FontSize',20)
ylabel('States','Interpreter','LaTeX','FontSize',20)
hold off

figure(2)
plot(t,x(:,8),'k-','LineWidth',2); 
hold on; 
plot(t,x(:,9),'r-','LineWidth',2);
legend({'$x_8=u_1$','$x_9=u_2$'},...
    'Interpreter','LaTeX','FontSize',18,'Location','NorthWest')
xlabel('Time (s)','Interpreter','LaTeX','FontSize',20)
ylabel('Controls','Interpreter','LaTeX','FontSize',20)
% subplot(2,2,1)
% plot(t,x(:,1),'o')
% xlabel('Time')
% ylabel('Horizontal Position')
% 
% subplot(2,2,2)
% plot(t,x(:,2),'o')
% xlabel('Time')
% ylabel('Vertical Position')
% 
% subplot(2,2,3)
% plot(t,x(:,3),'o')
% xlabel('Time')
% ylabel('Velocity')
% 
% subplot(2,2,4)
% plot(t,u,'o')
% xlabel('Time')
% ylabel('Control')