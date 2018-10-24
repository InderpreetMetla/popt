function output = lqrIHDynamics(t,x,u,auxdata)
tau = t;
x1 = x(:,1);
x2 = x(:,2); 
dt_dtau = 2./(1-tau).^2;

dx1 = dt_dtau.*x2;
dx2 = dt_dtau.*(2*x1-x2+u);
output = [dx1, dx2];
end