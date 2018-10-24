%******************************************************************%
% Free Flying Robot
% Sakawa, Y. (1999). "Trajectory Planning of a Free-Flying Robot 
% by Using the Optimal Control," Optimal Control Applications and  
% Methods", Vol. 20, pp. 235-248.
%
% It is also presented in Betts as Ex 6.13
%******************************************************************%
clear;clc
auxdata.alpha = 0.2;
auxdata.beta  = 0.2;
auxdata.gamma = 1;

t0 = 0;
tf = 12;
x0       = -10;  xf       = 0;
y0       = -10;  yf       = 0;
theta0   = pi/2; thetaf   = 0;
vx0      = 0;    vxf      = 0;
vy0      = 0;    vyf      = 0;
omega0   = 0;    omegaf   = 0;

xmin     = -10; xmax     = 10;
ymin     = -10; ymax     = 10;
thetamin = -pi; thetamax = pi;
vxmin    = -2;  vxmax    = 2;
vymin    = -2;  vymax    = 2;
omegamin = -1;  omegamax = 1;

%******************************************************************%
%  Setup for Problem Bounds  and Problem Guess
%******************************************************************%
u1Min = 0; u1Max = 1000;
u2Min = 0; u2Max = 1000;
u3Min = 0; u3Max = 1000;
u4Min = 0; u4Max = 1000;

bounds.phase.initialtime.lb = t0;
bounds.phase.initialtime.ub = t0;
bounds.phase.finaltime.lb   = tf;
bounds.phase.finaltime.ub   = tf;

bounds.phase.initialstate.lb = [x0,y0,theta0,vx0,vy0,omega0];
bounds.phase.initialstate.ub = [x0,y0,theta0,vx0,vy0,omega0];
bounds.phase.state.lb        = [xmin,ymin,thetamin,vxmin,vymin,omegamin];
bounds.phase.state.ub        = [xmax,ymax,thetamax,vxmax,vymax,omegamax];
bounds.phase.finalstate.lb   = [xf,yf,thetaf,vxf,vyf,omegaf];
bounds.phase.finalstate.ub   = [xf,yf,thetaf,vxf,vyf,omegaf];
bounds.phase.control.lb = [u1Min,u2Min,u3Min,u4Min];
bounds.phase.control.ub = [u1Max,u2Max,u3Max,u4Max];
bounds.phase.path.lb = [-1000, -1000];
bounds.phase.path.ub = [1, 1];

%******************************************************************%
% Setup Guess  
%******************************************************************%

tGuess     = [t0; tf];
xGuess     = [x0; xf];
yGuess     = [y0; yf];
thetaGuess = [theta0; thetaf];
vxGuess    = [vx0; vxf];
vyGuess    = [vy0; vyf];
omegaGuess = [omega0; omegaf];
u1Guess    = [0; 0];
u2Guess    = [0; 0];
u3Guess    = [0; 0];
u4Guess    = [0; 0];
guess.phase.time = tGuess;
guess.phase.state = [xGuess,yGuess,thetaGuess,vxGuess,vyGuess,omegaGuess];
guess.phase.control = [u1Guess,u2Guess,u3Guess,u4Guess];

%******************************************************************%
% Grid Information  
%******************************************************************%
grid.tol = 1e-6;
grid.phase.nodes.initialgrid = 2*ones(1,5);
grid.phase.nodes.lb = 9;
grid.phase.nodes.ub = 12;
grid.max_refine = 25;
%******************************************************************%
% Assemble Information into Problem Structure  
%******************************************************************%
problem.name                = 'FFR';
problem.funcs.Dynamics      = @FFRDynamics;
problem.funcs.PathObj       = @FFRPathObj;
problem.funcs.PathCst       = @FFRPathCst;
problem.funcs.BndObj        = @FFRBndObj; 
problem.funcs.BndCst        = @FFRBndCst;
problem.bounds              = bounds;
problem.guess               = guess;
problem.auxdata             = auxdata;
problem.derivatives.method	= 'ad';
problem.derivatives.order	= 1;
problem.grid                = grid;
problem.autoscale           = 'on';

%******************************************************************%
% Solve Problem
%******************************************************************%
[solution,plotdata]=popt(problem);
solution.cpu_time.total
%%
%******************************************************************%
% Plotting
%******************************************************************%
close all
t = solution.phase.time;
y = solution.phase.state;
u = solution.phase.control;

% Plot Solution
subplot(2,3,1)
plot(t,y(:,1),'k-','LineWidth',2.5);
xlabel('Time $(s)$','Interpreter','LaTeX','Fontsize',15);
ylabel('$x_1(t)~(m)$','Interpreter','LaTeX','Fontsize',15);
xlim([0 12])

subplot(2,3,2)
plot(t,y(:,2),'k-','LineWidth',2.5);
xlabel('Time $(s)$','Interpreter','LaTeX','Fontsize',15);
ylabel('$x_2(t)~(m)$','Interpreter','LaTeX','Fontsize',15);
xlim([0 12])

subplot(2,3,3)
plot(t,y(:,3),'k-','LineWidth',2.5);
xlabel('Time $(s)$','Interpreter','LaTeX','Fontsize',15);
ylabel('$x_3(t)~(rad)$','Interpreter','LaTeX','Fontsize',15);
xlim([0 12])

subplot(2,3,4)
plot(t,y(:,4),'k-','LineWidth',2.5);
xlabel('Time $(s)$','Interpreter','LaTeX','Fontsize',15);
ylabel('$x_4(t)~(m/s)$','Interpreter','LaTeX','Fontsize',15);
xlim([0 12])

subplot(2,3,5)
plot(t,y(:,5),'k-','LineWidth',2.5);
xlabel('Time $(s)$','Interpreter','LaTeX','Fontsize',15);
ylabel('$x_5(t)~(m/s)$','Interpreter','LaTeX','Fontsize',15);
xlim([0 12])

subplot(2,3,6)
plot(t,y(:,6),'k-','LineWidth',2.5);
xlabel('Time $(s)$','Interpreter','LaTeX','Fontsize',15);
ylabel('$x_6(t)~(rad/s)$','Interpreter','LaTeX','Fontsize',15);
xlim([0 12])

figure(2);
subplot(2,3,1)
plot(t,u(:,1),'ko-','LineWidth',1);
xlabel('Time $(s)$','Interpreter','LaTeX','Fontsize',15);
ylabel('$u_1(t)$','Interpreter','LaTeX','Fontsize',15);
xlim([0 12])

subplot(2,3,2)
plot(t,u(:,2),'k-','LineWidth',2.5);
xlabel('Time $(s)$','Interpreter','LaTeX','Fontsize',15);
ylabel('$u_2(t)$','Interpreter','LaTeX','Fontsize',15);
xlim([0 12])

subplot(2,3,3)
plot(t,u(:,3),'k-','LineWidth',2.5);
xlabel('Time $(s)$','Interpreter','LaTeX','Fontsize',15);
ylabel('$u_3(t)$','Interpreter','LaTeX','Fontsize',15);
xlim([0 12])

subplot(2,3,4)
plot(t,u(:,4),'k-','LineWidth',2.5);
xlabel('Time $(s)$','Interpreter','LaTeX','Fontsize',15);
ylabel('$u_4(t)$','Interpreter','LaTeX','Fontsize',15);
xlim([0 12])

subplot(2,3,5)
plot(t,u(:,1)-u(:,2),'k-','LineWidth',2.5);
xlabel('Time $(s)$','Interpreter','LaTeX','Fontsize',15);
ylabel('$T_1(t)$','Interpreter','LaTeX','Fontsize',15);
xlim([0 12])

subplot(2,3,6)
plot(t,u(:,3)-u(:,4),'k-','LineWidth',2.5);
xlabel('Time $(s)$','Interpreter','LaTeX','Fontsize',15);
ylabel('$T_2(t)$','Interpreter','LaTeX','Fontsize',15);
xlim([0 12])


for i=1:length(solution.grid_analysis.iteration)
    plot(solution.grid_analysis.iteration(i).phase.node_locations.',i,'ko')
    hold on
end