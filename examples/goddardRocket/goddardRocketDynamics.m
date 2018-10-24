function output = goddardRocketDynamics(t,x,u,auxdata)
h = x(:,1);
v = x(:,2);
m = x(:,3);
D = auxdata.dragk.*(v.^2).*exp(-h/auxdata.H);
hdot = v;
vdot = (u-D)./m-auxdata.g0*ones(size(t));
mdot = -u./auxdata.c;
output  = [hdot, vdot, mdot];
end