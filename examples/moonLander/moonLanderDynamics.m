function output = moonLanderDynamics(t,x,u,auxdata)
g = auxdata.g;
v = x(:,2);
dh = v;
dv = -g + u;
output = [dh, dv];
end
