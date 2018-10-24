function output = lqrIHPathObj(t,x,u,auxdata)
tau = t;
x1 = x(:,1);
x2 = x(:,2); 
dt_dtau = 2./(1-tau).^2;

output = dt_dtau.*(x1.^2+0.5*x2.^2+0.25*u.^2);
end