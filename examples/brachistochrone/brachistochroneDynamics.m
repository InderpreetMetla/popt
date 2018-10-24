function output = brachistochroneDynamics(t,x,u,auxdata)
g = auxdata.g;
v = x(:,3);
xdot = v.*sin(u);
ydot = v.*cos(u);
vdot = g*cos(u);
output = [xdot, ydot, vdot];
end