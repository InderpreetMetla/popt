%******************************************************************%
% Goddard Rocket Problem
% This problem is taken from the following reference:
% http://www.mcs.anl.gov/~more/cops/bcops/rocket.html
%******************************************************************%

clear;clc

h0 = 0;
v0 = 0;
m0 = 3;
mf = 1;
hmin = 0;
hmax = 30000;
vmin = -1500;
vmax = 1500;
mmin = 0.2*m0;
mmax = m0;
t0 = 0;
tfMin = 0;
tfMax = 500;

auxdata.g0 = 32.174;
auxdata.rho0 = 0.002378;
auxdata.H = 23800;
auxdata.csqrd = 3.264*auxdata.g0*auxdata.H;
auxdata.c = sqrt(auxdata.csqrd);
Tmax = 2*m0*auxdata.g0;
auxdata.dragk = 0.7110*Tmax/auxdata.csqrd;

%******************************************************************%
%  Setup for Problem Bounds  and Problem Guess
%******************************************************************%
%*********%
% Phase 1
%*********%
iphase = 1;
bounds.phase(iphase).initialtime.lb = t0;
bounds.phase(iphase).initialtime.ub = t0;
bounds.phase(iphase).finaltime.lb = tfMin;
bounds.phase(iphase).finaltime.ub = tfMax;
bounds.phase(iphase).initialstate.lb = [h0, v0, m0];
bounds.phase(iphase).initialstate.ub = [h0, v0, m0];
bounds.phase(iphase).state.lb = [hmin, vmin, mmin];
bounds.phase(iphase).state.ub = [hmax, vmax, mmax];
bounds.phase(iphase).finalstate.lb = [hmin, vmin, mmin];
bounds.phase(iphase).finalstate.ub = [hmax, vmax, mmax];
bounds.phase(iphase).control.lb = Tmax;
bounds.phase(iphase).control.ub = Tmax;
bounds.phase(iphase).path.lb = [];
bounds.phase(iphase).path.ub = [];
bounds.phase(iphase).boundary.lb = 0.1;
bounds.phase(iphase).boundary.ub = 1000;
% bounds.phase(iphase).nodes = 12; 
guess.phase(iphase).time =  [t0; tfMax];
guess.phase(iphase).state(:,1) = [h0; h0];
guess.phase(iphase).state(:,2) = [v0; v0];
guess.phase(iphase).state(:,3) = [m0; mf];
guess.phase(iphase).control = [0; Tmax];


%*********%
% Phase 2
%*********%
iphase = 2;
bounds.phase(iphase).initialtime.lb = tfMin;
bounds.phase(iphase).initialtime.ub = tfMax;
bounds.phase(iphase).finaltime.lb = tfMin;
bounds.phase(iphase).finaltime.ub = tfMax;
bounds.phase(iphase).initialstate.lb = [hmin, vmin, mmin];
bounds.phase(iphase).initialstate.ub = [hmax, vmax, mmax];
bounds.phase(iphase).state.lb = [hmin, vmin, mmin];
bounds.phase(iphase).state.ub = [hmax, vmax, mmax];
bounds.phase(iphase).finalstate.lb = [hmin, vmin, mmin];
bounds.phase(iphase).finalstate.ub = [hmax, vmax, mmax];
bounds.phase(iphase).control.lb = 0;
bounds.phase(iphase).control.ub = Tmax;
bounds.phase(iphase).path.lb = 0;
bounds.phase(iphase).path.ub = 0;
bounds.phase(iphase).boundary.lb = [0,0.1];
bounds.phase(iphase).boundary.ub = [0,1000];
% bounds.phase(iphase).nodes = 12; 
guess.phase(iphase).time =  [tfMin; tfMax];
guess.phase(iphase).state(:,1) = [h0; hmax];
guess.phase(iphase).state(:,2) = [vmax; vmax];
guess.phase(iphase).state(:,3) = [m0; mf];
guess.phase(iphase).control(:,1) = [0; Tmax];

%*********%
% Phase 3
%*********%
iphase = 3;
bounds.phase(iphase).initialtime.lb = tfMin;
bounds.phase(iphase).initialtime.ub = tfMax;
bounds.phase(iphase).finaltime.lb = tfMin;
bounds.phase(iphase).finaltime.ub = tfMax;
bounds.phase(iphase).initialstate.lb = [hmin, vmin, mmin];
bounds.phase(iphase).initialstate.ub = [hmax, vmax, mmax];
bounds.phase(iphase).state.lb = [hmin, vmin, mmin];
bounds.phase(iphase).state.ub = [hmax, vmax, mmax];
bounds.phase(iphase).finalstate.lb = [hmin, vmin, mf];
bounds.phase(iphase).finalstate.ub = [hmax, vmax, mf];
bounds.phase(iphase).control.lb = 0;
bounds.phase(iphase).control.ub = 0;
bounds.phase(iphase).path.lb = [];
bounds.phase(iphase).path.ub = [];
bounds.phase(iphase).boundary.lb = 0.1;
bounds.phase(iphase).boundary.ub = 1000;
% bounds.phase(iphase).nodes = 12; 
guess.phase(iphase).time =  [tfMin; tfMax];
guess.phase(iphase).state(:,1) = [hmin; hmax];
guess.phase(iphase).state(:,2) = [vmax; vmax];
guess.phase(iphase).state(:,3) = [m0; mf];
guess.phase(iphase).control(:,1) = [0; Tmax];

%******************************************************************%
%  Linkage Bounds
%******************************************************************%
ilink = 1;
bounds.link(ilink).lb = [0,0,0,0];
bounds.link(ilink).ub = [0,0,0,0];

ilink = 2;
bounds.link(ilink).lb = [0,0,0,0];
bounds.link(ilink).ub = [0,0,0,0];

%******************************************************************%
% Setup Grid
%******************************************************************%
grid.phase(1).nodes.initialgrid = 12*ones(1,2);
grid.phase(2).nodes.initialgrid = 12*ones(1,2);
grid.phase(3).nodes.initialgrid = 12*ones(1,2);
% grid.phase(1).nodes.initialgrid = 4*ones(1,1);
% grid.phase(2).nodes.initialgrid = 4*ones(1,1);
% grid.phase(3).nodes.initialgrid = 4*ones(1,1);
grid.tol = 1e-10;

%******************************************************************%
% Assemble Information into Problem Structure  
%******************************************************************%
problem.name                = 'goddardRocket';
problem.funcs.Dynamics      = @goddardRocketDynamics;
problem.funcs.PathObj       = @goddardRocketPathObj;
problem.funcs.PathCst       = @goddardRocketPathCst;
problem.funcs.BndObj        = @goddardRocketBndObj; 
problem.funcs.BndCst        = @goddardRocketBndCst;
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
tic
[solution, plotdata]   = popt(problem);
toc

%******************************************************************%
% Plotting
%******************************************************************%
%%
for iphase = 1:3
    t = plotdata.phase(iphase).time;
    x = plotdata.phase(iphase).state ;
    u = plotdata.phase(iphase).control;
    
    subplot(2,2,1)
    plot(t,x(:,1)/ 1000,'-','LineWidth',4)
    hold on
    xlabel('Time $(s)$','Interpreter','LaTeX','FontSize',15)
    ylabel('Height $(ft/1000)$','Interpreter','LaTeX','FontSize',15)
    
    subplot(2,2,2)
    plot(t,x(:,2),'-','LineWidth',4)
    hold on
    xlabel('Time $(s)$','Interpreter','LaTeX','FontSize',15)
    ylabel('Velocity $(ft/s)$','Interpreter','LaTeX','FontSize',15)
    
    subplot(2,2,3)
    plot(t,x(:,3),'-','LineWidth',4)
    hold on
    xlabel('Time $(s)$','Interpreter','LaTeX','FontSize',15)
    ylabel('Mass $(lb)$','Interpreter','LaTeX','FontSize',15)
    
    subplot(2,2,4)
    plot(t,u,'-','LineWidth',4)
    hold on
    xlabel('Time $(s)$','Interpreter','LaTeX','FontSize',15)
    ylabel('Control - Thrust $(lbf)$','Interpreter','LaTeX','FontSize',15)
    legend({'Phase 1', 'Phase 2', 'Phase 3'}, 'Interpreter','LaTeX','FontSize',12)
end