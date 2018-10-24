function output = lqrFHDynamics(t,x,u,auxdata)
xdot = 2*x(:,1) + 2*u(:,1).*(x(:,1).^(0.5));
output = xdot;
end