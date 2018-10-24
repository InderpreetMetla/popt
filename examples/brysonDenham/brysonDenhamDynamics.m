function output = brysonDenhamDynamics(t,x,u,auxdata)
v       = x(:,2);
xdot	= v;
vdot 	= u;
output  = [xdot, vdot];
end