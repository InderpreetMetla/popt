%----------------------- Moon-Lander Problem -----------------------------%
% This example can be found in the following reference:                   %
% Meditch, J., "On the Problem of Optimal Thrust Programming for a Soft   %
% Lunar Landing," IEEE Transactions on Automatic Control, Vol. 9,% No. 4, %
% 1964, pp. 477-484.                                                      %
%-------------------------------------------------------------------------%
clear all; clc

auxdata.g = 1.5;

t0min = 0;  t0max = 0;
tfmin = 0;  tfmax = 200;
h0 = 10;    hf = 0;
v0 = -2;    vf = 0;
hmin = 0;   hmax =  20;
vmin = -10; vmax =  10;
umin = 0;   umax = 3;

%******************************************************************%
%  Setup for Problem Bounds
%******************************************************************%
bounds.phase.initialtime.lb = t0min;
bounds.phase.initialtime.ub = t0max;
bounds.phase.finaltime.lb = tfmin;
bounds.phase.finaltime.ub = tfmax;
bounds.phase.initialstate.lb = [h0, v0];
bounds.phase.initialstate.ub = [h0, v0];
bounds.phase.state.lb = [hmin, vmin];
bounds.phase.state.ub = [hmax, vmax];
bounds.phase.finalstate.lb = [hf, vf];
bounds.phase.finalstate.ub = [hf, vf];
bounds.phase.control.lb = [umin];
bounds.phase.control.ub = [umax];

%******************************************************************%
% Provide Problem Guess
%******************************************************************% 
tGuess               = [t0min; 5];
hGuess               = [h0; hf];
vGuess               = [v0; vf];
uGuess               = [umin; umin];
guess.phase.state    = [hGuess, vGuess];
guess.phase.control  = [uGuess];
guess.phase.time     = [tGuess];

%******************************************************************%
% Setup Grid
%******************************************************************%
grid.phase.nodes.initialgrid = 8;
% grid.phase.nodes.initialgrid = [5,12,5,5];
grid.phase.nodes.lb = 4;
grid.phase.nodes.ub = 10;
grid.max_refine = 15;
grid.tol = 1e-6;

%******************************************************************%
% Assemble Information into Problem Structure  
%******************************************************************%
problem.name         	= 'Moon_Lander';
problem.funcs.Dynamics	= @moonLanderDynamics; 
problem.funcs.PathObj	= @moonLanderPathObj; 
problem.funcs.BndObj	= @moonLanderBndObj; 
problem.funcs.PathCst	= @moonLanderPathCst;
problem.funcs.BndCst	= @moonLanderBndCst;
problem.auxdata         = auxdata;
problem.bounds          = bounds;
problem.guess           = guess;
problem.grid            = grid;
problem.derivatives.method = 'fd';
problem.derivatives.order = 1;
problem.autoscale        = 'off';

%******************************************************************%
% Solve Problem
%******************************************************************%

solution = popt(problem);
solution.cpu_time.total

%******************************************************************%
% Plotting
%******************************************************************%
t = solution.phase.time;
x = solution.phase.state;
u = solution.phase.control;

% subplot(2,2,1)
% plot(t,x(:,1),'k-o')
% xlabel('Time')
% ylabel('Height')
% % grid on
% 
% subplot(2,2,2)
% plot(t,x(:,2),'k-o')
% xlabel('Time')
% ylabel('Velocity')
% % grid on
% 
% subplot(2,2,3)
% plot(t,u,'k-o')
% xlabel('Time')
% ylabel('Control')
% % grid on
% clear grid
subplot(1,2,1)
plot(t,x(:,1),'k-o','LineWidth',1.5)
hold on
plot(t,x(:,2),'k-v','LineWidth',1.5)
xlabel('Time (s)', 'Fontsize', 15)
ylabel('State', 'Fontsize', 15)
legend({'Height (m)', 'Velocity (m/s)'},'FontSize',12)
grid on
subplot(1,2,2)
plot(t,u(:,1),'k-o','LineWidth',1.5)
xlabel('Time (s)', 'Fontsize', 15)
ylabel('Control', 'Fontsize', 15)
legend({'Mass Flow Rate (kg/s)'},'FontSize',12)
grid on